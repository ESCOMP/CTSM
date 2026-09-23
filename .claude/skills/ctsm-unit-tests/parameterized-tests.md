# Running one CTSM test over several inputs

Companion to `SKILL.md`. Use when the same behaviour has to hold across several values of
a module-level configuration variable and you do not want the test body three times.

When the same behaviour has to hold across several values of a module-level configuration
variable, use pFUnit's parameterized test case. There is no example of this anywhere in
CTSM, so the whole pattern is below.

The alternatives are worse, and specifically so. A loop inside a single test body reports
**one** result no matter how many values it covers, and stops at the first failure —
pFUnit's assert macros expand to `...; if (anyExceptions()) return`, so the first bad
value returns out of the whole subroutine and the later ones never run. Copy-pasted
subroutines give you one result each but drift apart as the test evolves.

```fortran
module test_layer_count

  use funit
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar  , only : nlevsno
  use MyMod       , only : ComputeSomething

  implicit none
  save

  ! 1. The parameter type. toString is deferred in AbstractTestParameter,
  !    so you must implement it; it labels the run in the output.
  @testParameter
  type, extends(AbstractTestParameter) :: LayerCase
     integer :: nlev
   contains
     procedure :: toString
  end type LayerCase

  ! 2. The test type extends ParameterizedTestCase, NOT TestCase, and names
  !    both a constructor and the parameter list.
  @testCase(constructor=newTestLayerCount, testParameters={getParameters()})
  type, extends(ParameterizedTestCase) :: TestLayerCount
     integer :: nlev
   contains
     procedure :: setUp
     procedure :: tearDown
  end type TestLayerCount

contains

  ! 3. The constructor copies the parameter into the fixture instance.
  function newTestLayerCount(testParameter) result(test)
    type(LayerCase), intent(in) :: testParameter
    type(TestLayerCount)        :: test

    test%nlev = testParameter%nlev
  end function newTestLayerCount

  ! 4. The parameter list.
  function getParameters() result(params)
    type(LayerCase), allocatable :: params(:)

    allocate(params(3))
    params(1)%nlev = 3
    params(2)%nlev = 5
    params(3)%nlev = 12
  end function getParameters

  function toString(this) result(string)
    class(LayerCase), intent(in) :: this
    character(:), allocatable :: string
    character(len=32) :: buf

    write(buf, '(i0)') this%nlev
    string = 'nlev=' // trim(buf)
  end function toString

  ! setUp is how the parameter reaches the code under test: put it wherever
  ! the routine actually reads it from.
  subroutine setUp(this)
    class(TestLayerCount), intent(inout) :: this

    nlevsno = this%nlev
  end subroutine setUp

  subroutine tearDown(this)
    class(TestLayerCount), intent(inout) :: this
  end subroutine tearDown

  @Test
  subroutine computes_one_value_per_layer(this)
    class(TestLayerCount), intent(inout) :: this
    real(r8), allocatable :: result(:)

    call ComputeSomething(result)
    @assertEqual(this%nlev, size(result))
  end subroutine computes_one_value_per_layer

end module test_layer_count
```

Add the file to `pfunit_sources` in the directory's `CMakeLists.txt` like any other test.
Nothing else about the harness changes.

Omitting any of the four numbered pieces is a compile error rather than a silent
degradation, so the compiler will tell you. Two things it will not tell you:

- **Build the parameter list by assigning components, not with a structure constructor.**
  pFUnit's own examples write `LayerCase(3)`, which omits the componentless parent's
  component. It compiles under current `ifx`, but the explicit form above is the portable
  one.
- **`toString` is what makes a failure readable.** Keep it short; it is emitted twice.

### What the output looks like

A green run prints one undifferentiated dot per (test, parameter) and no names. Names
appear only with `-v`, and the `toString` label is emitted twice:

```
test_layer_count_suite.computes_one_value_per_layer[nlev=5][nlev=5]
```

A failure names the specific value that failed, gives the location, and **the remaining
parameter values still run**:

```
Failure in: test_layer_count_suite.computes_one_value_per_layer[nlev=5][nlev=5]
  Location: [test_layer_count.pf:86]
Tests run: 19, Failures: 1, Errors: 0
```

That last property is the whole reason to prefer this over a loop.
