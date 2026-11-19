---
name: designing-unit-test-cases
description: Use when choosing the values a test will set up, when every array position or both sides of an operation would get the same value, when a test's setup sits at a limit or is already consistent with the operation under test, when naming a test or writing its opening comment, when adding a sibling to a test that already exists, or when a mutation that should have been caught passes.
---

# Designing unit test cases

Two questions, asked of every test: **what class of bug is this setup blind to**, and **can a
reader who has never seen the code tell what the test pins down**.

**Core principle:** coverage comes from the *variation* in what a test sets up, not from the
number of assertions over it. A quantity that is the same everywhere pins nothing along the
dimension it is the same in, and no amount of asserting on it will change that.

This skill is about the values and the wording. The evidence a test owes — a red run or a
mutation — is in the `writing-tests-before-the-implementer` skill; the mechanics of building
and running are in `ctsm-unit-tests`.

**One word used throughout.** A **slot** is one position in an array the routine indexes into —
a soil or snow layer, a PFT, an aerosol species, a column — whether or not anything occupies it.
Where a CTSM routine appears below as an example, what that routine does is stated alongside it;
nothing here depends on already knowing it.

## Part 1 — Setups that do not test what you think

### The diagnostic question

For every dimension the setup spans — slot, layer, index, species, column, the two sides of
an operation, one call and the next — ask:

> **If the code were wrong along this dimension, would any number in the output change?**

If the answer is no, the test is blind there. A setup that holds a dimension constant, equal
or symmetric has traded that dimension's coverage away, usually without anyone noticing,
because the test still passes and still looks thorough.

Ask it **before** looking at the implementation. Never justify a setup by what the code does
today: these tests exist to constrain changes nobody has written yet, so "no current code path
does that" is a reason to expect the blind spot to matter later, not a reason to accept it.

### Five shapes that go blind

| Setup shape | What it cannot see | Instead |
|---|---|---|
| The same value in every slot | Any indexing error — reading the wrong slot, writing the wrong slot, or copying the array out and back with an offset all yield identical output | One distinct marker per slot |
| A two-index marker built as a product | A transposition of the two indices: `f(i)*g(j)` is symmetric, so swapped indices land on a value the setup already expects | A positional encoding, one field per index |
| A scale factor of one — or an offset of zero, an empty prefix, a weight of one | Both legs of a scale-and-unscale: at the identity element, applying and removing are each a no-op, so dropping either leg changes nothing | Any value that is not the operation's identity; a power of two keeps the arithmetic exact |
| Every case already at the limit | The limit: a setup that starts at the maximum behaves the same whether the limit is enforced, off by one, or absent | At least one case that has to grow *into* the limit, and one that stays below it |
| A setup the operation under test reproduces unchanged | That the operation stops where it should — running it one step too far rewrites the same values | A setup made inconsistent on purpose, so an overrun has something to overwrite |

### Markers, so a value landing in the wrong place is visible

One dimension. Distinct per slot, and never zero, so a marker cannot be mistaken for a slot the
routine zeroed:

```fortran
! j is the slot being filled and first_slot the lowest index of the array, so the
! multiplier runs 2, 3, 4, ... up the array. It starts at two rather than zero or one
! so that no marker collides with a zeroed slot or with an untouched count.
val = scale * real(j - first_slot + 2, r8)
```

Two dimensions — say a quantity carried per snow layer *and* per aerosol species, which the
routine reindexes one species at a time. **Do not build the marker as a product** — a product is
symmetric in its indices, so a reindex that transposed layer against species would land on
markers the setup already expects and pass unnoticed. Encode the two indices in separate
fields of one number, so every pair is distinct in both dimensions at once:

```fortran
! Species occupies the units digit, the slot multiplier everything above it. Only the
! species field is bounded, so the property survives a larger slot count: widening the
! slot field never reaches the units digit.
integer, parameter :: species_base = 10  ! must exceed the largest species index
val = scale * real(species_base * (j - first_slot + 2) + species, r8)
```

Give each quantity its own scale as well, so a value landing in the wrong *array* is visible
too, and write the resulting table into the comment — a reader should not have to evaluate the
expression to know what slot 3 holds.

### Limits: prefer growth to saturation

A setup that starts at a limit cannot demonstrate that the limit is enforced, and neither can
a set of cases that all sit below it.

CTSM's cap on how many snow layers a column may carry, `nlevsno`, is the concrete case. A
column set up with `nlevsno` layers already present never tries to exceed the cap, so it behaves
identically whether the cap is enforced, off by one, or missing altogether — and a test over it
passes in all three worlds. What pins the cap is a pair: one case that grows *into* it and is
held there, one that stays below it and is left alone.

The same holds for a floor, a clamp, and a maximum iteration count.

### The subtlest one: a setup that hides an overrun

Start with the concrete case. CTSM's vertical geometry is redundant: the interface at the top of
a layer is the interface at its bottom minus that layer's thickness, and the node sits at the
layer's midpoint. Routines that rearrange snow layers rebuild the nodes and interfaces from the
thicknesses afterwards, walking up the column one layer at a time, and stop when they reach the
top of the snow:

```fortran
do j = 0, -nlevsno+1, -1
   if (j >= snl(c)+1) then                    ! <- the stopping condition
      z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
      zi(c,j-1) = zi(c,j) - dz(c,j)
   end if
end do
```

Now suppose your test sets the column up so that its geometry is consistent everywhere —
including in the slots the loop is supposed to leave alone, where the node is also at its
layer's midpoint and the interfaces also agree with the thicknesses. Delete the stopping
condition and the loop runs one slot further. It computes a node depth and an interface depth
for that slot — and writes back exactly the numbers that were already there, because they were
consistent to begin with. Every assertion still passes, whether the routine stopped where it
should or overran by a slot — and so does a mutation that removes the stopping condition
outright. The coverage is zero, and it looks like coverage.

That is the general shape, and it has a name: the setup is a **fixed point** of the operation
under test when running that operation one step past where it should stop reproduces the setup
exactly. It matters wherever the point of a test is that something **stops** — a walk that must
halt at the edge of its region, a normalisation that must not run twice, a loop that must leave
the last element alone.

**Ask what one more step would change.** If the answer is nothing, then it is the *consistency*
of the setup that is hiding the overrun, and no assertion over a consistent setup can find it.

**So make the setup inconsistent on purpose.** Give the slot the routine must not touch a value
the operation would have produced differently — a node placed a quarter of the way down its layer
rather than at its midpoint, in the case above — so that an overrun has something to
overwrite and the assertion has something to catch. Then, at that value, say that it is
deliberately unphysical, what it buys, and why it cannot leak into anything else the test pins.
The safest place for it is a quantity nothing under test reads back — one the code only ever
writes.

### Holding a dimension constant on purpose

No test varies everything, and a test built to isolate one axis has to hold the others still.
That is a legitimate choice — but what makes it legitimate is that the blind spot it creates is
covered *somewhere in the file*, not that it is covered in this test.

So when you hold a dimension constant, do two things. Say in the detailed description that it is
deliberate and what it gives up, so the next reader does not take it for an oversight. And check
that some other test in the file varies it. A dimension that is constant in every test in the
file is not an isolation choice; it is exactly the blind spot the rest of Part 1 is about, with
a justification written over it.

## Part 2 — Tests a human can read

### Every test opens with a summary, then the detail

Three parts, in this order:

1. **One sentence.** It may wrap, and it may take a colon or a semicolon. It does not become
   two.
2. **A blank comment line.**
3. **The detailed description** — the numbers, the thresholds, why the setup holds the values
   it holds, and any cross-references.

The sentence has three properties:

- **It says what the routine under test must do in this case, and what makes the case worth
  pinning** — not what the test code does. `DivideSnowLayers` subdivides any snow layer that has
  grown thicker than the model allows a layer at its depth to be; so a test of the case where no
  layer has says "`DivideSnowLayers` must be the identity when no layer is over its maximum
  thickness", and never "sets up three layers and calls `DivideSnowLayers`". Name the
  configuration the case exercises: that is usually the axis the file is organised around, and
  it is not a number.
- **It contains no numbers.** No thicknesses, thresholds, layer indices or array bounds. In that
  same example, write "thicker than a bottom layer is allowed to be" rather than the thickness
  and the limit it exceeds. Numbers belong in the detailed description.
- **It stands alone.** No reference to another test, including implicit ones — not "the same
  setup as above", not "the complement of the previous test". A reader who jumps straight to
  one test gets the whole point from its first sentence. The detailed description below may
  cross-reference other tests freely; only the summary may not.

### Name a test for the condition and the branch, never for an input value

```
test_<routine>_<condition the setup creates>_<branch that condition drives>[_<axis>]
```

- **The condition is the state, not the number that produced it.** A number in a name means
  nothing to a reader who has not memorised the threshold it is meant to exceed, and it goes
  stale the moment that threshold moves.
- **The branch is what that condition drives** through the routine. A name that carries only
  the condition stops distinguishing anything the moment sibling tests exist for the other
  branches — and they will.
- **Keep the condition even when the branch name seems to cover it.** Which state the setup
  is in is often the whole reason one constraint binds rather than another; the branch name
  alone loses that.
- **Where a set of tests differs along one axis, name that axis in every one of them.** Not one
  test named for the subject and its twin named for the configuration: both carry the same
  suffix pair, so the difference is legible at a glance.

`DivideSnowLayers` again, which treats an over-thick layer differently depending on whether
anything lies beneath it — the bottom layer is halved, while a layer with another below it sheds
its excess downward. Two branches, so two conditions and two branch names:

```
test_divideSnowLayers_bottomLayerOverThickness_halved          ! condition + branch
test_divideSnowLayers_middleLayerOverThickness_shedsDownward   ! condition + branch
test_divideSnowLayers_dz0p2                                    ! names the input value
test_divideSnowLayers_geometry                                 ! names neither
```

**Expect to rename.** A name that distinguished a test while it stood alone stops distinguishing
once its siblings land, and renaming it then is the right move rather than a churn. Find out what
your test runner needs in order to notice the new name — in CTSM it needs a full rebuild, and the
check is that the count did *not* move; see `ctsm-unit-tests`.

### Give the quantities names, not literals

A bare `-4`, `-(my_nlevsno-1)` or `dz_layers(j+5)` forces the reader to re-derive the number's
meaning at every appearance, and the same value often means two different things in one test.
Declare named constants near the top and use them throughout — the number of occupied slots the
case starts and ends with, the first of them, any limit the case turns on — then write derived
indices in terms of them, so one expression replaces a different magic offset in every test.

Two payoffs beyond legibility:

- **The declaration is the natural place to say why the value is what it is** — including why a
  deliberately awkward value was chosen, which is where a fixed-point fix gets explained.
- **The differences between related tests collapse into the one line where the names are
  defined**, which is the whole point of the pair in miniature.

Where one value legitimately plays two roles, declare both names and let the second define the
first — `n_layers_at_start = n_layers_max`, in a case whose whole point is that the column
begins holding as many layers as it is allowed to hold. That documents the identity instead of
hiding it.

### Say what an array's dimensions are

Every `allocatable`, `pointer` or assumed-shape declaration in a test gets a trailing comment
naming its intended bounds:

```fortran
real(r8), allocatable :: snow_depth(:)     ! [begc:endc]
integer,  allocatable :: filter_snowc(:)   ! [1:num_snowc]
```

A bare `(:)` says nothing about which index space the test is working in — and where the test's
subject is indexing, the index space *is* the subject.

### Complementary tests announce each other and differ in as few lines as possible

Where a pair exists to isolate one axis, each **detailed description** names its twin and states
what is held fixed: "the same depth and the same layer count, differing only in whether the limit
is reached". Never in the topline, which has to stand alone.

Then hold the pair to it. A difference in setup order, in whether an assignment is spelled
`col%zi(c,0)` or `col%zi(bounds%begc:bounds%endc,0)`, or in the order of assertions is noise the
reader must rule out before finding the real difference. **Diff the two tests against each other
before committing them.**

## Red flags

- Every slot in the setup holds the same number
- A two-index marker built by multiplying the indices
- A fraction, weight or scale factor set to one, in a routine that both applies and removes it
- Every case sitting at the limit — or every case below it
- A mutation that should have been caught passes, and the conclusion is that the code is right
- A setup the routine under test would rewrite with the same numbers it already holds
- "No current code path does that" offered as a reason to hold a dimension constant
- A number, or an input value, in a test name
- Two sibling tests whose names disagree about which axis separates them
- A topline summary that carries a threshold value, or that starts "Tests that…" or "Sets up…"
- A topline that only makes sense after reading the test above it
- Two tests meant to differ in one thing, differing in five
- A bare `(:)` in a test declaration
