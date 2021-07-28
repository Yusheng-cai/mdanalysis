#pragma once
#include "XdrWrapper.h"

#include <string>
#include <iostream>

class TrrFile:public XdrWrapper
{
    public:
        TrrFile();
        virtual ~TrrFile(){};

        virtual void readNextFrame() override;
        virtual void readNumAtoms() override;    
};