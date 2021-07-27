#pragma once
#include "XdrWrapper.h"

#include <string>

class TrrFile:public XdrWrapper
{
    public:
        TrrFile():XdrWrapper(){};

        virtual void readNextFrame();

        virtual void readNumAtoms();    
};