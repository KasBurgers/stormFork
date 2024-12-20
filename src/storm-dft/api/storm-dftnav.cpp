#include "storm-dft/api/storm-dft.h"

#include "storm-conv/api/storm-conv.h"
#include "storm-conv/settings/modules/JaniExportSettings.h"
#include "storm-dft/settings/modules/DftGspnSettings.h"
#include "storm-dft/settings/modules/FaultTreeSettings.h"

#include <memory>
#include <vector>
#include "storm-dft/adapters/SFTBDDPropertyFormulaAdapter.h"
#include "storm-dft/modelchecker/DftModularizationChecker.h"
#include "storm-dft/modelchecker/SFTBDDChecker.h"
#include "storm-dft/storage/DFT.h"
#include "storm-dft/storage/DftJsonExporter.h"
#include "storm-dft/storage/SylvanBddManager.h"
#include "storm-dft/transformations/SftToBddTransformator.h"
#include "storm-dft/utility/MTTFHelper.h"

namespace storm::dft {
namespace api {

template<>
void analyzeDFTBdd(std::shared_ptr<storm::dft::storage::DFT<double>> const& dft, bool const exportToDot, std::string const& filename, bool const calculateMttf,
                   double const mttfPrecision, double const mttfStepsize, std::string const mttfAlgorithmName, bool const calculateMCS,
                   bool const calculateProbability, bool const useModularisation, std::string const importanceMeasureName,
                   std::vector<double> const& timepoints, std::vector<std::shared_ptr<storm::logic::Formula const>> const& properties,
                   std::vector<std::string> const& additionalRelevantEventNames, size_t const chunksize) {
    if (calculateMttf) {
        if (mttfAlgorithmName == "proceeding") {
            std::cout << "The numerically approximated MTTF is " << storm::dft::utility::MTTFHelperProceeding(dft, mttfStepsize, mttfPrecision) << '\n';
        } else if (mttfAlgorithmName == "variableChange") {
            std::cout << "The numerically approximated MTTF is " << storm::dft::utility::MTTFHelperVariableChange(dft, mttfStepsize) << '\n';
        }
    }

    if (useModularisation && calculateProbability) {
        storm::dft::modelchecker::DftModularizationChecker checker{dft};
        if (chunksize == 1) {
            for (auto const& timebound : timepoints) {
                auto const probability{checker.getProbabilityAtTimebound(timebound)};
                std::cout << "System failure probability at timebound " << timebound << " is " << probability << '\n';
            }
        } else {
            auto const probabilities{checker.getProbabilitiesAtTimepoints(timepoints, chunksize)};
            for (size_t i{0}; i < timepoints.size(); ++i) {
                auto const timebound{timepoints[i]};
                auto const probability{probabilities[i]};
                std::cout << "System failure probability at timebound " << timebound << " is " << probability << '\n';
            }
        }
        if (!properties.empty()) {
            auto const probabilities{checker.check(properties, chunksize)};
            for (size_t i{0}; i < probabilities.size(); ++i) {
                std::cout << "Property \"" << properties.at(i)->toString() << "\" has result " << probabilities.at(i) << '\n';
            }
        }
        return;
    } else {
        STORM_LOG_THROW(dft->nrDynamicElements() == 0, storm::exceptions::NotSupportedException,
                        "DFT is dynamic. "
                        "Bdds can only be used on static fault trees. "
                        "Try modularisation.");
    }

    auto sylvanBddManager{std::make_shared<storm::dft::storage::SylvanBddManager>()};
    sylvanBddManager->execute([&]() {
        storm::dft::utility::RelevantEvents relevantEvents{additionalRelevantEventNames.begin(), additionalRelevantEventNames.end()};
        storm::dft::adapters::SFTBDDPropertyFormulaAdapter adapter{dft, properties, relevantEvents, sylvanBddManager};
        auto checker{adapter.getSFTBDDChecker()};

        if (exportToDot) {
            checker->exportBddToDot(filename);
        }

        if (calculateMCS) {
            auto const minimalCutSets{checker->getMinimalCutSetsAsIndices()};
            auto const sylvanBddManager{checker->getSylvanBddManager()};

            std::cout << "{\n";
            for (auto const& minimalCutSet : minimalCutSets) {
                std::cout << '{';
                for (auto const& be : minimalCutSet) {
                    std::cout << sylvanBddManager->getName(be) << ' ';
                }
                std::cout << "},\n";
            }
            std::cout << "}\n";
        }

        if (calculateProbability & (!dft->getOrderedValuePairOfSwitchGates().empty()) & 0) {
            std::cout << "Start family-based analysis " << '\n';
            auto listOfPairs = dft->getOrderedValuePairOfSwitchGates();
            auto subModules = dft->getAllSubModules();
            auto remainingFT = *dft;
            auto remainingFTOut = *dft;
            std::list<storm::dft::storage::DFT<double>> ftList;
            std::deque<std::string> orderedListOfNames;
            std::deque<double> output;

            for (auto x : listOfPairs) {
                auto id = x.first;
                if (dft->isRootOfSubtree(id)) {
                    auto name = remainingFT.getElement(id)->name();
                    orderedListOfNames.push_back(name);
                }
            }

            auto topLevelName = remainingFT.getElement(remainingFT.getTopLevelIndex())->name();
            // make sure topLevel is analysed even when top is no switch
            if (std::find(orderedListOfNames.begin(), orderedListOfNames.end(), topLevelName) == orderedListOfNames.end()) {
                orderedListOfNames.push_back(topLevelName);
            }
            ftList.push_back(remainingFTOut);
            orderedListOfNames.push_back("placeholder");

            while (!orderedListOfNames.empty()) {
                // if (ftList.empty()) {
                orderedListOfNames.pop_back();
                //}
                // update remainingFT

                auto ftListCopy = ftList;
                for (auto differentFT : ftListCopy) {
                    remainingFT = differentFT;

                    // update modules
                    subModules = remainingFT.getAllSubModules();
                    bool swapped = false;

                    for (auto subFt : subModules) {
                        // get name of representative
                        auto ID = subFt.getRepresentative();
                        auto elementName = remainingFT.getElement(ID)->name();
                        // break loop early if list is empty already
                        if (orderedListOfNames.empty()) {
                            break;
                        }

                        // debug purpuse
                        // std::cout << "value 1 =  " << orderedListOfNames.back() << '\n';
                        // std::cout << "value 2 =  " << elementName << '\n';

                        if (orderedListOfNames.back() == elementName && swapped == false) {
                            // std::cout << "inside ifloop 2 " << '\n';
                            // pop current element
                            if (ftList.empty()) {
                                orderedListOfNames.pop_back();
                            } else {
                                ftList.pop_front();
                            }
                            // swap has occured no need to swap again this itteration
                            swapped = true;

                            auto newDft = subFt.getSubtree(remainingFT);
                            storm::dft::adapters::SFTBDDPropertyFormulaAdapter adapter2{std::make_shared<storm::dft::storage::DFT<double>>(newDft), properties,
                                                                                        relevantEvents, sylvanBddManager};
                            auto checker2{adapter2.getSFTBDDChecker()};
                            std::string outputProb = "";
                            std::list<double> probList;
                            for (auto const& timebound : timepoints) {
                                auto const probability{checker2->getProbabilityAtTimeboundSwitch(timebound)};
                                for (auto x : probability) {
                                    outputProb = std::to_string(x);
                                    // std::cout << outputProb << '\n';
                                }
                                probList = probability;
                            }
                            for (auto x : probList) {
                                outputProb = std::to_string(x);
                                remainingFTOut = newDft.replaceSubtree(remainingFT, subFt, outputProb);

                                if (!remainingFTOut.isBasicElement(remainingFTOut.getTopLevelIndex())) {
                                    ftList.push_back(remainingFTOut);
                                } else {
                                    output.push_back(x);
                                }
                            }
                        }
                        if (swapped) {
                            break;
                        }
                    }
                    if (!swapped) {
                        orderedListOfNames.pop_back();
                    }
                }
            }
            for (auto out : output) {
                std::cout << "Prob " << out << '\n';
            }
            // make MTBDD out of leaves
            std::deque<storm::dd::Add<storm::dd::DdType::Sylvan, double>> terminal_nodes;
            std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
            std::deque<std::string> MTBDDorderedListOfNames;

            auto startFT = *dft;
            for (auto x : listOfPairs) {
                auto id = x.first;
                auto name = startFT.getElement(id)->name();
                MTBDDorderedListOfNames.push_back(name);
            }

            std::cout << "Order: ";
            for (auto name : MTBDDorderedListOfNames) {
                std::cout << name << ", ";
            }
            std::cout << '\n';

            for (int i = 0; i < output.size(); i++) {
                terminal_nodes.push_back(manager->getConstant(output[i]));
            }
            std::deque<std::pair<storm::expressions::Variable, storm::expressions::Variable>> variableQueue;
            for (int i = 0; i < MTBDDorderedListOfNames.size(); i++) {
                variableQueue.push_back(manager->addMetaVariable(MTBDDorderedListOfNames[i], 0, 1));
            }

            std::deque<storm::dd::Add<storm::dd::DdType::Sylvan, double>> nodeQueue;
            for (auto variable : variableQueue) {
                // make nodes containing terminal nodes
                if (variable == variableQueue[0]) {
                    for (int i = 0; i < output.size(); i += 2) {
                        nodeQueue.push_back(manager->getEncoding(variable.first, 0).ite(terminal_nodes[i], terminal_nodes[i + 1]));
                    }
                }
                // first nodes are made
                else {
                    auto loopSize = nodeQueue.size() / 2;
                    for (int i = 0; i < loopSize; i++) {
                        nodeQueue.push_back(manager->getEncoding(variable.first, 0).ite(nodeQueue[0], nodeQueue[1]));
                        // pop used nodes
                        nodeQueue.pop_front();
                        nodeQueue.pop_front();
                    }
                }
            }
            // nodeQueue[0].exportToDot(filename + "MTBB.txt");
        }

        // if (calculateProbability & (dft->getOrderedValuePairOfSwitchGates().empty())) {
        if (calculateProbability) {
            std::cout << "Start naive family-based analysis " << '\n';
            storm::utility::Stopwatch NaiveTimer;
            NaiveTimer.start();

            if (chunksize == 1) {
                for (auto const& timebound : timepoints) {
                    auto const probability{checker->getProbabilityAtTimebound(timebound)};
                    std::cout << "System failure probability at timebound " << timebound << " is " << probability << '\n';
                }
            } else {
                auto const probabilities{checker->getProbabilitiesAtTimepoints(timepoints, chunksize)};
                for (size_t i{0}; i < timepoints.size(); ++i) {
                    auto const timebound{timepoints[i]};
                    auto const probability{probabilities[i]};
                    std::cout << "System failure probability at timebound " << timebound << " is " << probability << '\n';
                }
            }

            if (!properties.empty()) {
                auto const probabilities{adapter.check(chunksize)};
                for (size_t i{0}; i < probabilities.size(); ++i) {
                    std::cout << "Property \"" << properties.at(i)->toString() << "\" has result " << probabilities.at(i) << '\n';
                }
            }

            NaiveTimer.stop();
            std::cout << "Naive time " << NaiveTimer.getTimeInMilliseconds() << "ms" << '\n';
        }

        if (importanceMeasureName != "" && timepoints.size() == 1) {
            auto const bes{dft->getBasicElements()};
            std::vector<double> values{};
            if (importanceMeasureName == "MIF") {
                values = checker->getAllBirnbaumFactorsAtTimebound(timepoints[0]);
            }
            if (importanceMeasureName == "CIF") {
                values = checker->getAllCIFsAtTimebound(timepoints[0]);
            }
            if (importanceMeasureName == "DIF") {
                values = checker->getAllDIFsAtTimebound(timepoints[0]);
            }
            if (importanceMeasureName == "RAW") {
                values = checker->getAllRAWsAtTimebound(timepoints[0]);
            }
            if (importanceMeasureName == "RRW") {
                values = checker->getAllRRWsAtTimebound(timepoints[0]);
            }

            for (size_t i{0}; i < bes.size(); ++i) {
                std::cout << importanceMeasureName << " for the basic event " << bes[i]->name() << " at timebound " << timepoints[0] << " is " << values[i]
                          << '\n';
            }
        } else if (importanceMeasureName != "") {
            auto const bes{dft->getBasicElements()};
            std::vector<std::vector<double>> values{};
            if (importanceMeasureName == "MIF") {
                values = checker->getAllBirnbaumFactorsAtTimepoints(timepoints, chunksize);
            }
            if (importanceMeasureName == "CIF") {
                values = checker->getAllCIFsAtTimepoints(timepoints, chunksize);
            }
            if (importanceMeasureName == "DIF") {
                values = checker->getAllDIFsAtTimepoints(timepoints, chunksize);
            }
            if (importanceMeasureName == "RAW") {
                values = checker->getAllRAWsAtTimepoints(timepoints, chunksize);
            }
            if (importanceMeasureName == "RRW") {
                values = checker->getAllRRWsAtTimepoints(timepoints, chunksize);
            }
            for (size_t i{0}; i < bes.size(); ++i) {
                for (size_t j{0}; j < timepoints.size(); ++j) {
                    std::cout << importanceMeasureName << " for the basic event " << bes[i]->name() << " at timebound " << timepoints[j] << " is "
                              << values[i][j] << '\n';
                }
            }
        }
    });
}

template<>
void analyzeDFTBdd(std::shared_ptr<storm::dft::storage::DFT<storm::RationalFunction>> const& dft, bool const exportToDot, std::string const& filename,
                   bool const calculateMttf, double const mttfPrecision, double const mttfStepsize, std::string const mttfAlgorithmName,
                   bool const calculateMCS, bool const calculateProbability, bool const useModularisation, std::string const importanceMeasureName,
                   std::vector<double> const& timepoints, std::vector<std::shared_ptr<storm::logic::Formula const>> const& properties,
                   std::vector<std::string> const& additionalRelevantEventNames, size_t const chunksize) {
    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "BDD analysis is not supportet for this data type.");
}

template<typename ValueType>
void exportDFTToJsonFile(storm::dft::storage::DFT<ValueType> const& dft, std::string const& file) {
    storm::dft::storage::DftJsonExporter<ValueType>::toFile(dft, file);
}

template<typename ValueType>
std::string exportDFTToJsonString(storm::dft::storage::DFT<ValueType> const& dft) {
    std::stringstream stream;
    storm::dft::storage::DftJsonExporter<ValueType>::toStream(dft, stream);
    return stream.str();
}

template<>
void exportDFTToSMT(storm::dft::storage::DFT<double> const& dft, std::string const& file) {
    storm::dft::modelchecker::DFTASFChecker asfChecker(dft);
    asfChecker.convert();
    asfChecker.toFile(file);
}

template<>
void exportDFTToSMT(storm::dft::storage::DFT<storm::RationalFunction> const& dft, std::string const& file) {
    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Export to SMT does not support this data type.");
}

template<>
void analyzeDFTSMT(storm::dft::storage::DFT<double> const& dft, bool printOutput) {
    uint64_t solverTimeout = 10;

    storm::dft::modelchecker::DFTASFChecker smtChecker(dft);
    smtChecker.toSolver();
    // Removed bound computation etc. here
    smtChecker.setSolverTimeout(solverTimeout);
    smtChecker.checkTleNeverFailed();
    smtChecker.unsetSolverTimeout();
}

template<>
void analyzeDFTSMT(storm::dft::storage::DFT<storm::RationalFunction> const& dft, bool printOutput) {
    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Analysis by SMT not supported for this data type.");
}

template<>
std::pair<std::shared_ptr<storm::gspn::GSPN>, uint64_t> transformToGSPN(storm::dft::storage::DFT<double> const& dft) {
    storm::dft::settings::modules::FaultTreeSettings const& ftSettings = storm::settings::getModule<storm::dft::settings::modules::FaultTreeSettings>();
    storm::dft::settings::modules::DftGspnSettings const& dftGspnSettings = storm::settings::getModule<storm::dft::settings::modules::DftGspnSettings>();

    // Set Don't Care elements
    std::set<uint64_t> dontCareElements;
    if (!ftSettings.isDisableDC()) {
        // Insert all elements as Don't Care elements
        for (std::size_t i = 0; i < dft.nrElements(); i++) {
            dontCareElements.insert(dft.getElement(i)->id());
        }
    }

    // Transform to GSPN
    storm::dft::transformations::DftToGspnTransformator<double> gspnTransformator(dft);
    auto priorities = gspnTransformator.computePriorities(dftGspnSettings.isExtendPriorities());
    gspnTransformator.transform(priorities, dontCareElements, !dftGspnSettings.isDisableSmartTransformation(), dftGspnSettings.isMergeDCFailed(),
                                dftGspnSettings.isExtendPriorities());
    std::shared_ptr<storm::gspn::GSPN> gspn(gspnTransformator.obtainGSPN());
    return std::make_pair(gspn, gspnTransformator.toplevelFailedPlaceId());
}

std::shared_ptr<storm::jani::Model> transformToJani(storm::gspn::GSPN const& gspn, uint64_t toplevelFailedPlace) {
    // Build Jani model
    storm::builder::JaniGSPNBuilder builder(gspn);
    std::shared_ptr<storm::jani::Model> model(builder.build("dft_gspn"));

    // Build properties
    std::shared_ptr<storm::expressions::ExpressionManager> const& exprManager = gspn.getExpressionManager();
    storm::jani::Variable const& topfailedVar = builder.getPlaceVariable(toplevelFailedPlace);
    storm::expressions::Expression targetExpression = exprManager->integer(1) == topfailedVar.getExpressionVariable().getExpression();
    // Add variable for easier access to 'failed' state
    builder.addTransientVariable(model.get(), "failed", targetExpression);
    auto failedFormula = std::make_shared<storm::logic::AtomicExpressionFormula>(targetExpression);
    auto properties = builder.getStandardProperties(model.get(), failedFormula, "Failed", "a failed state", true);

    // Export Jani to file
    storm::dft::settings::modules::DftGspnSettings const& dftGspnSettings = storm::settings::getModule<storm::dft::settings::modules::DftGspnSettings>();
    if (dftGspnSettings.isWriteToJaniSet()) {
        auto const& jani = storm::settings::getModule<storm::settings::modules::JaniExportSettings>();
        storm::api::exportJaniToFile(*model, properties, dftGspnSettings.getWriteToJaniFilename(), jani.isCompactJsonSet());
    }

    return model;
}

template<>
std::pair<std::shared_ptr<storm::gspn::GSPN>, uint64_t> transformToGSPN(storm::dft::storage::DFT<storm::RationalFunction> const& dft) {
    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Transformation to GSPN not supported for this data type.");
}

// Explicitly instantiate methods
template void exportDFTToJsonFile(storm::dft::storage::DFT<double> const&, std::string const&);
template std::string exportDFTToJsonString(storm::dft::storage::DFT<double> const&);

template void exportDFTToJsonFile(storm::dft::storage::DFT<storm::RationalFunction> const&, std::string const&);
template std::string exportDFTToJsonString(storm::dft::storage::DFT<storm::RationalFunction> const&);

}  // namespace api
}  // namespace storm::dft
