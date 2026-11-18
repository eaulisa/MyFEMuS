These runners reuse all 50 C++ sources and their original inputs/assertions from `unittests/testCutFemWeightLSNative`, with the original class and the local ex47 class in separate executables.

Build both backends with the same compiler and optimization:

```sh
python applications/MultiPhase/ex47/native_tests/build_comparison.py --build /tmp/ex47-native --boost /path/to/boost --eigen /path/to/eigen3 --jobs 2
python applications/MultiPhase/ex47/native_tests/build_comparison.py --build /tmp/ex47-native --boost /path/to/boost --eigen /path/to/eigen3 --jobs 2 --tests native_suite_speed curved_public_regression
```

Run public-class accuracy checks and all other fixtures:

```sh
python applications/MultiPhase/ex47/native_tests/run_native_cases.py --build /tmp/ex47-native --output /tmp/ex47-native-results/native-final --jobs 2
python applications/MultiPhase/ex47/native_tests/run_comparison.py --build /tmp/ex47-native --output /tmp/ex47-native-results/accuracy --group all --jobs 2
```

After compilation and other runs finish, collect interleaved timing batches:

```sh
python applications/MultiPhase/ex47/native_tests/time_native_cases.py --build /tmp/ex47-native --accuracy /tmp/ex47-native-results/native-final --output /tmp/ex47-native-results/native-timing --samples 3
python applications/MultiPhase/ex47/native_tests/run_comparison.py --build /tmp/ex47-native --output /tmp/ex47-native-results/fixture-timing --group focused --jobs 1 --samples 3
```

`MeasuredWeights.hpp` measures construction separately from weight calls and records quadrature work. The speed driver clears final weights before each repeated public API call; its original accuracy assessor still checks every timed result. These batches recompute quadrature while retaining one rule and its prepared element/geometry data. The existing ex47 benchmarks additionally alternate two retained element records and exercise all surface modes, including WEDGE21.

`OriginalCompatibility.cpp` forwards primitive calls to production. `NativeCompatibility.cpp` translates their original native/unit coordinate conventions and routes integration into the prepared kernels; it never samples polynomial inputs back into double nodal values. Timings of these adapters include compatibility overhead and are distinct from public-class timing. Coefficient import is limited to layouts supported by the data class. Utility methods inherited by compatibility classes (power integrals, chart algebra, and other unchanged helpers) remain shared infrastructure. `planar_suite` exercises the legacy planar class and is an auxiliary reference, not a second native implementation.

The original detector uses depth 12. The fixture runner records the new default depth 1 separately and repeats detection with `EX47_NATIVE_DETECTION_DEPTH=12`; conservative uncertified cells are not interpreted as incorrect cut integrals.

The additional `curved_public_regression` reuses the perturbed curved HEX27 fixture at degree 2 through the actual public weight APIs. The unit-coordinate correction resolves that failure; both backends are retested through the same public API. Known failing fixtures intentionally retain their original assertions and produce a nonzero runner status.

The timing-only graph and degree-13 surface drivers also have accuracy-only variants that expose every moment without adding work to the timed executables:

```sh
python applications/MultiPhase/ex47/native_tests/build_comparison.py --build /tmp/ex47-native --boost /path/to/boost --eigen /path/to/eigen3 --jobs 2 --tests benchmark_graph_accuracy benchmark_surface_accuracy
python applications/MultiPhase/ex47/native_tests/run_benchmark_accuracy.py --build /tmp/ex47-native --output /tmp/ex47-native-results/benchmark-accuracy
python applications/MultiPhase/ex47/native_tests/run_ctest_comparison.py --build /tmp/ex47-native --output /tmp/ex47-native-results/ctest --detection-depth 12
python applications/MultiPhase/ex47/native_tests/make_report.py --build /tmp/ex47-native --results /tmp/ex47-native-results --output applications/MultiPhase/ex47/native_tests/results/unit_coordinates --unit-coordinates
```

The polynomial-coefficient import supports the primitive fixtures without converting their long-double inputs into nodal doubles. Public prepared-weight overloads continue to require a root record created from FEM node values.

Unit-coordinate validation is recorded separately under `results/unit_coordinates/`; the other result tables are the historical baseline before this correction.
