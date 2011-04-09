#include <string.h>

#include "vorticity.h"

extern "C" int test_op(int argc, char ** argv) {
	SphereOperator op(19, 36, 0);
	SphereVorticity vor(op);

	bool test_vor = false;
	bool res = true;

	if (argc == 0) {
		test_vor = true;
	} else {
		for (int i = 0; i < argc; ++i) {
			if (!strcmp(argv[i], "vorticity")) {
				test_vor = true;
			}
		}
	}

	if (test_vor) {
		res &= vor.test();
	}

	if (res) {
		return 0;
	} else {
		return -1;
	}
}
