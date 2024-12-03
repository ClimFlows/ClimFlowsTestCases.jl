var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ClimFlowsTestCases","category":"page"},{"location":"#ClimFlowsTestCases","page":"Home","title":"ClimFlowsTestCases","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ClimFlowsTestCases.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ClimFlowsTestCases]","category":"page"},{"location":"#ClimFlowsTestCases.DCMIP","page":"Home","title":"ClimFlowsTestCases.DCMIP","text":"struct DCMIP{N} <: TestCaseHPE\n\nDCMIP 2012 test cases. Currently only N=42 (test case 4.2) is implemented.\n\n\n\n\n\n","category":"type"},{"location":"#ClimFlowsTestCases.Isothermal","page":"Home","title":"ClimFlowsTestCases.Isothermal","text":"struct Isothermal <: TestCaseHPE\n\nReturn testcase corresponding to an initial state of rest with exponential decrease of pressure with geopotential. For a perfect gas, this coincides with the classic isothermal profile, and psurf * vsurf = R*T with T the constant temperature. The surface geopotential defines a wiggly mountain. Numerics based on terrain-following coordinates typically struggle with such case and generate spurious motion.\n\n\n\n\n\n","category":"type"},{"location":"#ClimFlowsTestCases.Jablonowski06","page":"Home","title":"ClimFlowsTestCases.Jablonowski06","text":"struct Jablonowski06{P} <: TestCaseHPE\n\nJablonowski (2006) and DCMIP 4.2 test cases (dry vs. moist). The only difference is specific humidity is included for the moist case.\n\n\n\n\n\n","category":"type"},{"location":"#ClimFlowsTestCases.TestCase","page":"Home","title":"ClimFlowsTestCases.TestCase","text":"abstract type TestCase end\n\nParent type for test case types. Derived types must implement a constructor accepting an optional precision and keyword arguments to override default ones. For example:\n\ncase = Williams91{6}(Float32; R0=3e-3)\n\nDerived types must also possess a field params which is a named tuple containing the parameters, and implement initial and describe. \n\n\n\n\n\n","category":"type"},{"location":"#ClimFlowsTestCases.TestCaseHPE","page":"Home","title":"ClimFlowsTestCases.TestCaseHPE","text":"Parent type for hydrostatic test cases.\n\n\n\n\n\n","category":"type"},{"location":"#ClimFlowsTestCases.TestCaseSW","page":"Home","title":"ClimFlowsTestCases.TestCaseSW","text":"Parent type for shallow-water test cases.\n\n\n\n\n\n","category":"type"},{"location":"#ClimFlowsTestCases.Williamson91","page":"Home","title":"ClimFlowsTestCases.Williamson91","text":"struct Williamson91{N} <: TestCaseSW\n\nWilliamson (1991) test cases. Currently only N=6 is implemented.\n\n\n\n\n\n","category":"type"},{"location":"#ClimFlowsTestCases.describe","page":"Home","title":"ClimFlowsTestCases.describe","text":"string = describe(case::TestCase)\n\nPlain text description of case, possibly multi-line.\n\n\n\n\n\n","category":"function"},{"location":"#ClimFlowsTestCases.initial","page":"Home","title":"ClimFlowsTestCases.initial","text":"gh, ulon, ulat = initial(case::TestCaseSW, lon, lat)\n\nFor a shallow-water test case, returns the initial geopotential thickness and velocity at given latitude and longitude.\n\nps, Phis = initial(case::TestCaseHPE, lon, lat)\n\nFor a hydrostatic test case, returns the initial surface pressure and geopotential at given latitude and longitude.\n\nPhi, ulon, ulat, q = initial(case::TestCaseHPE, lon, lat, p)\n\nFor a hydrostatic test case, returns the initial geopotential, velocity and composition at given latitude, longitude and pressure.\n\n\n\n\n\n","category":"function"}]
}
