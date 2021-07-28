#include "XdrWrapper.h"

class XtcFile:public XdrWrapper
{
    public:
        XtcFile();
        virtual ~XtcFile(){};

        virtual void readNextFrame() override;
        virtual void readNumAtoms() override;
    private:
        Frame::Real precision_;
};