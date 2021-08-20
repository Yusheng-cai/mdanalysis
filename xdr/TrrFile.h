#pragma once
#include "XdrWrapper.h"
#include "libxdr/trr_seek.h"

#include <string>
#include <iostream>

class TrrFile:public XdrWrapper
{
    public:
        TrrFile(const XdrInput& input);
        virtual ~TrrFile(){};

        virtual bool readFrame(int FrameNum) override;
        virtual void readNumAtoms() override;    
        virtual void readNframes() override;
    private:
        std::vector<int> offsets_;
};