#pragma once
#include "XdrWrapper.h"

#include <string>
#include <iostream>

class TrrFile:public XdrWrapper
{
    public:
        TrrFile(const XdrInput& input);
        virtual ~TrrFile(){};

        virtual bool readNextFrame() override;
        virtual void readNumAtoms() override;    
};