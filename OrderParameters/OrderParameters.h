#include "tools/GenericFactory.h"
#include "tools/InputParser.h"

class OrderParameters
{
    public:
        OrderParameters();
        ~OrderParameters();

        virtual void calculate();
    private:
};