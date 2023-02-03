#include "Calculation.h"
#include "SRE.h"
#include "Bin.h"

#include <memory>
#include <vector>

class SlabVsr : public Calculation{
    using SRE_ptr = std::unique_ptr<SRE>;

    public:
        SlabVsr(const CalculationInput& input);

        virtual void calculate() override;

        virtual void update() override{};

    private:
        SRE_ptr vsr_;
};