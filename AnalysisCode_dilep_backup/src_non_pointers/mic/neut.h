
#include <vector>
#include <cmath>
#include <sys/time.h>
#include "../myvector.h"
#include "../dilep_input.h"
#include "../utilities.h"
#include "calc.h"

#define TO1D(nc,tid,sol,did) nc[tid*16+sol*4+did]

namespace Dilep {
	namespace Phi {
		void dilep (std::vector<DilepInput> &vdi);
	}
}
