#pragma once
#include "XdrWrapper.h"

#include <string>
#include <iostream>

class TrrFile:public XdrWrapper
{
    public:
        TrrFile();
        virtual ~TrrFile(){};

        virtual bool readNextFrame() override;
        virtual void readNumAtoms() override;    
};