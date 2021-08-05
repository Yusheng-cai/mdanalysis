#include "XdrWrapper.h"
#include "libxdr/xtc_seek.h"

class XtcFile:public XdrWrapper
{
    public:
        XtcFile(const XdrInput& input);
        virtual ~XtcFile(){};

        virtual bool readNextFrame() override;
        virtual void readNumAtoms() override;
        virtual void readNframes() override;
    private:
        Frame::Real precision_;
};