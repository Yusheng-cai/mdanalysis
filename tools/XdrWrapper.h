#include "libxdr/xdrfile_trr.h"
#include "libxdr/xdrfile_xtc.h"
#include "libxdr/xdrfile.h"

#include <string>

class XdrWrapper
{
    public:
        XdrWrapper();    
        ~XdrWrapper();
    private:
        XDRFILE* file;
        int natoms_;
};