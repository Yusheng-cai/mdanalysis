#pragma once

#include "tools/InputParser.h"
#include "tools/CommonTypes.h"

#include <vector>
#include <functional>
#include <string>

class OutputValue
{
    public:
        using Real          = CommonTypes::Real;

        using VectorReal    = std::vector<Real>;
        using ValueFunction = std::function<Real(void)>;

        OutputValue(std::string name, const ValueFunction& func):name_(name),val_func_(func){};

        Real get_value() const {return val_func_();};
        std::string get_name() const {return name_;};

    private:
        std::string name_;
        const ValueFunction val_func_;
};