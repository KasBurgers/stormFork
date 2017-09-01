#include "storm/modelchecker/multiobjective/pcaa/SparseMdpPcaaWeightVectorChecker.h"

#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/utility/macros.h"
#include "storm/utility/vector.h"
#include "storm/logic/Formulas.h"
#include "storm/solver/MinMaxLinearEquationSolver.h"
#include "storm/solver/LinearEquationSolver.h"


#include "storm/exceptions/InvalidPropertyException.h"
#include "storm/exceptions/IllegalArgumentException.h"
#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"


namespace storm {
    namespace modelchecker {
        namespace multiobjective {
            
            template <class SparseMdpModelType>
            SparseMdpPcaaWeightVectorChecker<SparseMdpModelType>::SparseMdpPcaaWeightVectorChecker(SparseMultiObjectivePreprocessorResult<SparseMdpModelType> const& preprocessorResult) :
                SparsePcaaWeightVectorChecker<SparseMdpModelType>(preprocessorResult) {
                this->initialize(preprocessorResult);
            }
            
            template <class SparseMdpModelType>
            void SparseMdpPcaaWeightVectorChecker<SparseMdpModelType>::initializeModelTypeSpecificData(SparseMdpModelType const& model) {
                // set the state action rewards. Also do some sanity checks on the objectives.
                this->actionRewards.resize(this->objectives.size());
                for (uint_fast64_t objIndex = 0; objIndex < this->objectives.size(); ++objIndex) {
                    auto const& formula = *this->objectives[objIndex].formula;
                    STORM_LOG_THROW(formula.isRewardOperatorFormula() && formula.asRewardOperatorFormula().hasRewardModelName(), storm::exceptions::UnexpectedException, "Unexpected type of operator formula: " << formula);
                    STORM_LOG_THROW(formula.getSubformula().isCumulativeRewardFormula() || formula.getSubformula().isTotalRewardFormula(), storm::exceptions::UnexpectedException, "Unexpected type of sub-formula: " << formula.getSubformula());
                    typename SparseMdpModelType::RewardModelType const& rewModel = model.getRewardModel(formula.asRewardOperatorFormula().getRewardModelName());
                    STORM_LOG_THROW(!rewModel.hasTransitionRewards(), storm::exceptions::NotSupportedException, "Reward model has transition rewards which is not expected.");
                    this->actionRewards[objIndex] = rewModel.getTotalRewardVector(model.getTransitionMatrix());
                }
            }
            
            template <class SparseMdpModelType>
            void SparseMdpPcaaWeightVectorChecker<SparseMdpModelType>::boundedPhase(std::vector<ValueType> const& weightVector, std::vector<ValueType>& weightedRewardVector) {
                // Allocate some memory so this does not need to happen for each time epoch
                std::vector<uint_fast64_t> optimalChoicesInCurrentEpoch(this->transitionMatrix.getRowGroupCount());
                std::vector<ValueType> choiceValues(weightedRewardVector.size());
                std::vector<ValueType> temporaryResult(this->transitionMatrix.getRowGroupCount());
                // Get for each occurring timeBound the indices of the objectives with that bound.
                std::map<uint_fast64_t, storm::storage::BitVector, std::greater<uint_fast64_t>> stepBounds;
                for (uint_fast64_t objIndex = 0; objIndex < this->objectives.size(); ++objIndex) {
                    if (this->objectives[objIndex].formula->getSubformula().isCumulativeRewardFormula()) {
                        auto const& subformula = this->objectives[objIndex].formula->getSubformula().asCumulativeRewardFormula();
                        uint_fast64_t stepBound = subformula.template getBound<uint_fast64_t>();
                        if (subformula.isBoundStrict()) {
                            --stepBound;
                        }
                        auto stepBoundIt = stepBounds.insert(std::make_pair(stepBound, storm::storage::BitVector(this->objectives.size(), false))).first;
                        stepBoundIt->second.set(objIndex);
                        
                        // There is no error for the values of these objectives.
                        this->offsetsToUnderApproximation[objIndex] = storm::utility::zero<ValueType>();
                        this->offsetsToOverApproximation[objIndex] = storm::utility::zero<ValueType>();
                    }
                }
                
                // Stores the objectives for which we need to compute values in the current time epoch.
                storm::storage::BitVector consideredObjectives = this->objectivesWithNoUpperTimeBound;
                
                auto stepBoundIt = stepBounds.begin();
                uint_fast64_t currentEpoch = stepBounds.empty() ? 0 : stepBoundIt->first;
                
                while (currentEpoch > 0) {
                    
                    if (stepBoundIt != stepBounds.end() && currentEpoch == stepBoundIt->first) {
                        consideredObjectives |= stepBoundIt->second;
                        for(auto objIndex : stepBoundIt->second) {
                            // This objective now plays a role in the weighted sum
                            ValueType factor = storm::solver::minimize(this->objectives[objIndex].formula->getOptimalityType()) ? -weightVector[objIndex] : weightVector[objIndex];
                            storm::utility::vector::addScaledVector(weightedRewardVector, this->actionRewards[objIndex], factor);
                        }
                        ++stepBoundIt;
                    }
                    
                    // Get values and scheduler for weighted sum of objectives
                    this->transitionMatrix.multiplyWithVector(this->weightedResult, choiceValues);
                    storm::utility::vector::addVectors(choiceValues, weightedRewardVector, choiceValues);
                    storm::utility::vector::reduceVectorMax(choiceValues, this->weightedResult, this->transitionMatrix.getRowGroupIndices(), &optimalChoicesInCurrentEpoch);
                    
                    // get values for individual objectives
                    // TODO we could compute the result for one of the objectives from the weighted result, the given weight vector, and the remaining objective results.
                    for (auto objIndex : consideredObjectives) {
                        std::vector<ValueType>& objectiveResult = this->objectiveResults[objIndex];
                        std::vector<ValueType> const& objectiveRewards = this->actionRewards[objIndex];
                        auto rowGroupIndexIt = this->transitionMatrix.getRowGroupIndices().begin();
                        auto optimalChoiceIt = optimalChoicesInCurrentEpoch.begin();
                        for(ValueType& stateValue : temporaryResult){
                            uint_fast64_t row = (*rowGroupIndexIt) + (*optimalChoiceIt);
                            ++rowGroupIndexIt;
                            ++optimalChoiceIt;
                            stateValue = objectiveRewards[row];
                            for(auto const& entry : this->transitionMatrix.getRow(row)) {
                                stateValue += entry.getValue() * objectiveResult[entry.getColumn()];
                            }
                        }
                        objectiveResult.swap(temporaryResult);
                    }
                    --currentEpoch;
                }
            }
            
            template class SparseMdpPcaaWeightVectorChecker<storm::models::sparse::Mdp<double>>;
#ifdef STORM_HAVE_CARL
            template class SparseMdpPcaaWeightVectorChecker<storm::models::sparse::Mdp<storm::RationalNumber>>;
#endif
            
        }
    }
}
