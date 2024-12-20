#pragma once

#include "DFTGate.h"

namespace storm::dft {
namespace storage {
namespace elements {

/*!
 * SWITCH gate.
 * Fails if at least one child has failed.
 */
template<typename ValueType>
class DFTSwitch : public DFTGate<ValueType> {
   public:
    /*!
     * Constructor.
     * @param id Id.
     * @param name Name.
     * @param children Children.
     * @param decisionVariable Decision variable used to choose left or right subtree.
     */
    DFTSwitch(size_t id, std::string const& name, std::vector<std::shared_ptr<DFTElement<ValueType>>> const& children = {})
        : DFTGate<ValueType>(id, name, children) {
    }

    std::shared_ptr<DFTElement<ValueType>> clone() const override {
        return std::shared_ptr<DFTElement<ValueType>>(new DFTSwitch<ValueType>(this->id(), this->name(), {}));
    }

    storm::dft::storage::elements::DFTElementType type() const override {
        return storm::dft::storage::elements::DFTElementType::SWITCH;
    }

    bool isStaticElement() const override {
        return false;
    }

    void checkFails(storm::dft::storage::DFTState<ValueType>& state, storm::dft::storage::DFTStateSpaceGenerationQueues<ValueType>& queues) const override {
    // TODO
    
    }

    void checkFailsafe(storm::dft::storage::DFTState<ValueType>& state, storm::dft::storage::DFTStateSpaceGenerationQueues<ValueType>& queues) const override {
    // TODO

    // All children have failed
            this->fail(state, queues);
    }
};

}  // namespace elements
}  // namespace storage
}  // namespace storm::dft
