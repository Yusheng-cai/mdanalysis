#pragma once
#include "XdrWrapper.h"

#include <string>
#include <iostream>

class TrrFile:public XdrWrapper
{
    public:
        TrrFile(const ParameterPack& pack);
        virtual ~TrrFile(){};

        virtual bool readNextFrame() override;
        virtual void readNumAtoms() override;    
};