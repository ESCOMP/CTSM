---
orphan: true
---

(clm-coding-conventions-page)=

# CLM coding conventions

(clm-coding-conventions)=

## CLM Coding Conventions

This page is intended as a guide for CLM coding standards, practices and conventions. On this page, we focus on five items that we want all developers to concentrate on. Keep these in the back of your mind as you are programming and refer to them if you have a question. We are more than willing to work with you throughout the process as needed; please contact us at [ctsm-support@ucar.edu](mailto:ctsm-support@ucar.edu) or on your GitHub discussion, issue, or pull request. If code doesn't generally conform to these items then it will delay the tagging process and we will ask you to modify your code accordingly.

Our past experiences have led us to document a more [Comprehensive list of standards](#comprehensive-list-of-standards) that we are striving for. You will find constructs that do not follow these guidelines and we are updating the legacy CLM code base as time and resources permit.

### Comments

1. Comment often and include explanations for algorithms and processes. Your code should be easily readable and understandable by anyone else who reads your work.

```{rubric} Logical Branches
:heading-level: 4
```

include a comment after an **else** or **endif** statement. Use something that makes sense for the type of logical branch you are using. This is particularly important when there are large numbers of nested blocks. Look at the second if statement below for an even nicer way to accomplish the same thing without comments.

```{code-block} fortran
:caption: Comments and Logical Branches

if ( someLogical ) then ! someLogical controls the assignment of x or y
   x = 0.0_r8
else ! else someLogical
   y = 0.0_r8
endif ! if someLogical

! Alternative to comments
fates_if: if (use_fates) then
   lch4_if: if (use_lch4) then
[...]
   end if lch4_if
end if fates_if
```

```{rubric} Pointer Intent
:heading-level: 4
```

If you pass information to a subroutine via a pointer in an associate block, include a comment as to the **intent** (OUT, IN, INOUT) of that pointer. This is in addition to a description of that pointer as well as the type and array dimensions. Also see [**Argument Passing**](#clm-conventions-argument-passing) on this page.

```{code-block} fortran
:caption: Comments and Pointer Intent

associate(&
   ptrIn    => derivedType%readOnlyInRoutine      , & ! [character(len=8) (:)] intent = in    used for x
   ptrInOut => derivedType%readAndWriteInRoutine  , & ! [real(r8) (:,:,:)]     intent = inout used for y
   ptrOut   => derivedType%writeOnlyInRoutine       & ! [real(r8) (:)]         intent = out   used for z
)
...
end associate
```

1. The meaning (semantic) of **intent** (out, in, inout)
   1. out - a variable is only set within a subroutine
   2. in - a variable is only passed into and used by a subroutine
   3. inout - a variable may be both modified and used within a subroutine

(clm-conventions-argument-passing)=

### Argument Passing

1. The **intent** attribute of subroutine input/output variables is used. Argument lists should ***always*** use the **intent** attribute (try to have data be either IN or OUT and avoid INOUT if possible/reasonable).
2. There are currently three ways in the current code with which you can pass information to subroutines. By interface (best choice), by pointer (second best choice) and by using a global data type (a bad option, but sometimes still necessary within the current CLM framework). We are in the process of refactoring the code so that all data will be passed through the subroutine interface. New code coming in should not pass via global data unless all other options do not work.

```{rubric} Subroutine Interface
:heading-level: 4
```

Use the intent keyword and keyword pairs if you have more than a few arguments. Try to consistently pass either all derived types or all FORTRAN types. Passing derived types can shorten the list of arguments and make it more understandable. Below are two examples from CLM, ([subsection i](#clm-conventions-clm-example), below) and RTM, ([subsection ii](#clm-conventions-rtm-example), below), respectively.

(clm-conventions-clm-example)=

```{rubric} CLM example from BandDiagonal::BandDiagonalMod.F90
:heading-level: 5
```

```{rubric} Subroutine definition
:heading-level: 6
```

```{code-block} text
:caption: Interfaces::CLM. Example from BandDiagonal, definition

subroutine BandDiagonal(bounds, lbj, ubj, jtop, jbot, numf, filter, nband, b, r, u)
    !
    ! !DESCRIPTION:
    ! Tridiagonal matrix solution
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                    ! bounds
    integer , intent(in)    :: lbj, ubj                        ! lbinning and ubing level indices
    integer , intent(in)    :: jtop( bounds%begc: )            ! top level for each column [col]
    integer , intent(in)    :: jbot( bounds%begc: )            ! bottom level for each column [col]
    integer , intent(in)    :: numf                            ! filter dimension
    integer , intent(in)    :: nband                           ! band width
    integer , intent(in)    :: filter(:)                       ! filter
    real(r8), intent(in)    :: b( bounds%begc: , 1:   , lbj: ) ! compact band matrix [col, nband, j]
    real(r8), intent(in)    :: r( bounds%begc: , lbj: )        ! "r" rhs of linear system [col, j]
    real(r8), intent(inout) :: u( bounds%begc: , lbj: )        ! solution [col, j]
    !
    ! ! LOCAL VARIABLES:
    integer  :: j,ci,fc,info,m,n              !indices
    integer  :: kl,ku                         !number of sub/super diagonals
    integer, allocatable :: ipiv(:)           !temporary
    real(r8),allocatable :: ab(:,:),temp(:,:) !compact storage array
    real(r8),allocatable :: result(:)

    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(jtop) == (/bounds%endc/)),             errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(jbot) == (/bounds%endc/)),             errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(b)    == (/bounds%endc, nband, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(r)    == (/bounds%endc, ubj/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(u)    == (/bounds%endc, ubj/)),        errMsg(__FILE__, __LINE__))
```

The definition above illustrates a few important points.
- Specify the lower bounds of array arguments (this is needed for the subroutine to operate properly in a threaded region)
- Do **not** specify the upper bounds of array arguments (doing so can prevent the compiler from checking for array size agreement between the caller and callee, and can carry a performance penalty)
- The first set of code in the routine should be assertions for the expected array upper bounds

The above code demonstrates a few stylistic points as well:
- Specify the lower bound explicitly, even when it is 1 (see the 'b' argument, above)
- Put a space after each ':' for readability
- Add comments about expected array sizes, as in \[col, nband, j\]

```{rubric} Subroutine call
:heading-level: 6
```

Follow this example in which BandDiagonal is called from SoilTemperatureMod.F90:

```{code-block} fortran
:caption: Interfaces::CLM. Example from BandDiagonal, call

call BandDiagonal(bounds, -nlevsno, nlevgrnd, jtop(begc:endc), jbot(begc:endc), &
         num_nolakec, filter_nolakec, nband, bmatrix(begc:endc, :, :), &
         rvector(begc:endc, :), tvector(begc:endc, :))
```

- The important point here is that the lower and upper bounds are explicitly specified for the gridcell / landunit / col / pft dimension of all array arguments. This is important for threading to work properly.
- Even for local variables whose dimensions only go begc:endc, please explicitly specify the dimensions as above when passing the array to a subroutine. This will ensure that threading continues to work if the variable's dimensions change (e.g., if it is changed to have dimensions that span the proc_bounds rather than the clump_bounds).

(clm-conventions-rtm-example)=

```{rubric} RTM example from RtmMod.F90
:heading-level: 5
```

```{code-block} text
:caption: Interfaces::RTM. Example from RtmMod.F90

!
! The call of the routine
!
call RtmFloodInit (frivinp_rtm, begr, endr, nt_rtm, &
     runoff%fthresh( begr:endr ),                   &
     evel( begr:endr , : ),                         &
     runoff%gindex( begr:endr ),                    &
     runoff%lnumr,                                  &
     flood_active,                                  &
     effvel_active)
...
!
! The interface definition
!
subroutine RtmFloodInit(frivinp, begr, endr, nt_rtm, fthresh, evel, &
     gindex ,        &
     lnumr ,         &
     is_rtmflood_on, &
     is_effvel_on    )
...
   ! Subroutine arguments
   ! in mode arguments
   character(len=*), intent(in) :: frivinp
   integer ,         intent(in) :: begr, endr, nt_rtm
   logical ,         intent(in) :: is_rtmflood_on      !control flooding
   logical ,         intent(in) :: is_effvel_on        !control eff. velocity
   integer ,         intent(in) :: gindex( begr: )     ! global index [begr:endr]
   integer ,         intent(in) :: lnumr               ! local number of cells
...
! check bounds of arrays
   SHR_ASSERT_ALL((ubound(fthresh) == (/endr/)),         errMsg(__FILE__, __LINE__))
   SHR_ASSERT_ALL((ubound(gindex)  == (/endr/)),         errMsg(__FILE__, __LINE__))
   SHR_ASSERT_ALL((ubound(evel)    == (/endr, nt_rtm/)), errMsg(__FILE__, __LINE__))
...
```

```{rubric} Pointer
:heading-level: 4
```

We are moving away from using pointers as a way pass information into routines. We still allow this but with the use of an associate statement.

```{code-block} fortran
:caption: Comments and Pointer Intent

associate(&
   ptrIn    => derivedType%readOnlyInRoutine     , & ! [character(len=8) (:)] intent = in used for x
   ptrInOut => derivedType%readAndWriteInRoutine , & ! [real(r8) (:,:,:)]     intent = inout used for y
   ptrOut   => derivedType%writeOnlyInRoutine      & ! [real(r8) (:)]         intent = out used for z
)
...
end associate
```

```{rubric} Global data
:heading-level: 4
```

We are moving away from using global data as a way pass information into routines. Unless you are extending existing functionality (e.g. using a namelist variable to control some execution) we will generally not allow any new use of a global variable. The other exception to this rule is that it is OK to you use global type declarations.

```{code-block} fortran
:caption: Global Data. Global data in CanopyFluxesMod.F90

use clm_varctl  , only: iulog, use_cn, use_lch4, use_c13, use_c14, use_cndv
!
   if ( use_c13 ) then
      c13o2(p) = forc_pc13o2(g)
   endif ! use_c13
!
```

```{code-block} fortran
:caption: Global Data. Using a global type declaration

use decompMod   , only : bounds_type ! using a global type declaration
...
subroutine CNSoilLittVertTransp(bounds, num_soilc, filter_soilc)
...
type(bounds_type), intent(in) :: bounds  ! bounds
...
```

### Indentation

1. Indent consistently and make sure your auto-tab function in your editor converts tabs to spaces.
   1. Indent of most blocks (do, if, while, etc.): 3 spaces
   2. Indent of "program" blocks (program, module, subroutine, function): 2 spaces
   3. Indent of "associate" blocks: 0 spaces
   4. Indent of continuation lines: 5 (this can vary based on look and feel)

   ```{code-block} fortran
   :caption: Indentation

   ! a dash, '-', represents a space
   --use thisModule , only : thisRoutine
   !
   --subroutine foo
   ----use thatModule , only : thatRoutine
   ----associate ( &
   -------ptr1 => someDerivedType%ptr1 , &
   -------ptr2 => someDerivedType%ptr2   &
   ----)
   ----if ( someLogical ) then
   -------x=x+1
   ! leave a line before and after statements that stand alone.
   -------if (x == 0) then endrun()
   ! leave a line before and after statements that stand alone.
   -------do while ( x <= 1000 )
   ----------! do some stuff
   ----------z = 1.0_r8 / radius
   -------end do ! end x<=1000
   ----endif ! someLogical
   !
   ----end associate
   --end subroutine
   ```
2. Do not leave hanging white space at the end of any line of code; git diff shows this in red spaces.
3. Do not tabularize code by adding extra spaces to any line of code. It may be tempting to do this for better readability, but it also takes time to manage code that needs realignment when new code comes in.

### Statements and Subroutines

1. Don't repeat yourself. If you are copying and pasting code, then write a subroutine or function to encapsulate functionality.

   ```{code-block} fortran
   :caption: Example of repeated code - Do

   !
   ! an example of using a subprogram to remove repeated code where
   ! the max daylength is calculated in a function.
   !
   ! From iniTimeConst:

         ! initialize maximum daylength, based on latitude and maximum declination
         ! maximum declination hardwired for present-day orbital parameters,
         ! ± 23.4667 degrees = ± 0.409571 radians, use negative value for S. Hem
         max_decl = 0.409571
         if (grc%lat(g) .lt. 0._r8) max_decl = -max_decl
            max_dayl(c) = calcMaxDayFunc( grc%lat(g) , max_decl )

   ! From CanopyFluxesMod:

         ! calculate daylength
         dayl(c) = calcMaxDayFunc( lat(g) , decl(c) )

   ! From CNPhenologyMod (copied & pasted in two places in this module):

            lat = (SHR_CONST_PI/180._r8)*grc%latdeg(pgridcell(p))
            dayl(p) = calcMaxDayFunc( lat , decl(c) )
   ```

   ```{code-block} fortran
   :caption: Example of repeated code - Don't do

   !
   ! an example of what not to do in terms of repeating code.
   ! The max daylength calculation is repeated.
   !
   ! From iniTimeConst:

         ! initialize maximum daylength, based on latitude and maximum declination
         ! maximum declination hardwired for present-day orbital parameters,
         ! ± 23.4667 degrees = ± 0.409571 radians, use negative value for S. Hem
         max_decl = 0.409571
         if (grc%lat(g) .lt. 0._r8) max_decl = -max_decl
            temp = -(sin(grc%lat(g))*sin(max_decl))/(cos(grc%lat(g)) * cos(max_decl))
            temp = min(1._r8,max(-1._r8,temp))
            max_dayl(c) = 2.0_r8 * 13750.9871_r8 * acos(temp)

   ! From CanopyFluxesMod:

         ! calculate daylength
         temp = -(sin(lat(g))*sin(decl(c)))/(cos(lat(g)) * cos(decl(c)))
         temp = min(1._r8,max(-1._r8,temp))
         dayl = 2.0_r8 * 13750.9871_r8 * acos(temp)

   ! From CNPhenologyMod (copied & pasted in two places in this module):

            lat = (SHR_CONST_PI/180._r8)*grc%latdeg(pgridcell(p))
            temp = -(sin(lat)*sin(decl(c)))/(cos(lat) * cos(decl(c)))
            temp = min(1._r8,max(-1._r8,temp))
            dayl(p) = 2.0_r8 * 13750.9871_r8 * acos(temp)
   ```
2. Fortran allows you to place multiple statements on one line. Don't do this.

   ```{code-block} fortran
   :caption: Multiple Statements - Do

   ! do this:
   x = 0
   y = 0
   z = 0
   ```

   ```{code-block} fortran
   :caption: Multiple Statements - Don't do

   ! not this:
   x = 0 ; y=0; z=0
   ```
3. Use temporary variables often

   ```{code-block} fortran
   :caption: Temporary Variables - Do

   ! do this:
   someFactor = ( const1/const2 ) * const3
   additiveFluxes = ( add1Flx + add2Flx + add3Flx ) ** someFactor
   subFluxes = (sub1Flx + sub2Flx + sub3Flx) ** someFactor
   sumVal = additiveFluxes + subFluxes
   ```

   ```{code-block} fortran
   :caption: Temporary Variables - Don't do

   ! not this:
   sumVal = (( add1Flx + add2Flx + add3Flx ) ** (( const1/const2 ) * const3 )) + ((sub1Flx + sub2Flx + sub3Flx) ** (( const1/const2 ) * const3 ))
   ```
4. Keep subroutines short and organize the code into the smallest reasonable chunk of code that can stand on its own. Basically, you're looking for dividing routines into logical units of work. For example, don't put initialization, some science and then diagnostics all in one routine, these would be best divided into 4 routines; in this case one routine that calls three additional routines.
   - Use parenthesis to clarify you statements even when the Fortran precedence rules will evaluate items in the correct order.

### Preprocessor Macros (CPP Tokens)

- We are in the process of removing CPP Tokens from all of CLM (see [our refactoring page](https://wiki.ucar.edu/spaces/ccsm/pages/274339715/Recent+CLM4.5+Refactoring#RecentCLM4.5Refactoring-Preprocessormacroremovalmacro)) and will not accept code using them.

  ```{code-block} fortran
  :caption: Preprocessor Macros (CPP Tokens) - Do

  ! do this:
  if (useNewPhysics == .true.) then
     ...
     ! do some new physics
     ...
  endif
  ```

  ```{code-block} fortran
  :caption: Preprocessor Macros (CPP Tokens) - Don't do

  ! not this:
  #if (defined USENEWPHYSICS)
     ...
     ! do some new physics
     ...
  #endif
  ```

(comprehensive-list-of-standards)=

## Comprehensive list of standards

These conventions cover some of the more detailed aspects of CLM code design. It is intended to compliment and extend the information found on the [previous page](#clm-coding-conventions). If you have questions please contact us at [ctsm-support@ucar.edu](mailto:ctsm-support@ucar.edu) or on your GitHub discussion, issue, or pull request.

### General Guidelines

1. We support [Fortran 2003](http://fortranwiki.org/fortran/show/Fortran+2003+status) language features with a few exceptions, along with some [Fortran 2008](http://fortranwiki.org/fortran/show/Fortran+2008+status) features. The Fortran 2003 features that we don't implement are parameterized derived types and derived type I/O. The Fortran 2008 features that we do implement are `newunit=` in `open`, the `erf` and `gamma` intrinsics, and `storage_size`.
2. Use **implicit none**. All variables must be explicitly typed.
3. **Don't repeat yourself**. Duplicated code is difficult to maintain. If you find yourself using copy and paste more than once, create a new function or subroutine that encapsulates the desired functionality.
4. Try to implement subroutine and function interface with argument lists of 9 or less. If you have lists longer than this, break up the functionality into multiple methods and/or use derived types for higher level interfaces. Try to have argument lists short enough that you can make unit-tests for them.
5. If you have argument lists greater than 8 -- you must call such routines using the keyword-pair form for all arguments (but see 3 above).

### Performance, debugging and optimization

1. Don't mix up C/C++ and Fortran memory access order; it will slow your code down drastically.
2. Debug with a debugger (totalview, gdb, idb, pgdbg, etc...) it will save you a lot of time finding errors.
   1. Adding Write statements to your code and dumping variable values also can work very well in finding problems.
3. Be aware of the memory footprint of additional code. Try to keep increases in memory use low. Use [Valgrind](http://valgrind.org/) or similar tools often during your development process to track and monitor memory increases, changes, leaks, etc....
4. Performance
   1. Regarding vectorization. If you feel your code performance is suffering there are some ideas [here](http://software.intel.com/en-us/intel-vectorization-tools) under the Fortran section.
   2. How you access derived type members can change performance. See for example some ideas regarding [derived types](http://software.intel.com/sites/products/documentation/doclib/stdxe/2013/composerxe/compiler/fortran-lin/index.htm#GUID-AFE9908C-BDC7-47DE-BEE6-2C7FBA382F4B.htm) and a discussion on [alignment](http://software.intel.com/sites/products/documentation/doclib/stdxe/2013/composerxe/compiler/fortran-mac/GUID-E9F21218-F280-4521-86B2-FAE99FF0BA0F.htm)

### Modularity

1. Data types, data and methods should be encapsulated within a module if they are used in that module. If data types, data and/or methods are used outside of the module they are declared in, they must be made publicly available.
2. Data should also be accessible at the most private level possible. In other words:
   1. If a variable is only used in one subroutine it should be local to that subroutine
   2. If a variable is only local to one module then it is declared private to that module.
   3. Public data for modules should be used sparingly.
3. Model state and flux variables live in derived types defined in their own modules (e.g., `CropType.F90`). Add new variables to the relevant type rather than to shared global modules.
4. Derived types provide their own `Init`, `InitHistory`, `InitCold`, and `Restart` procedures as necessary.

### Guidelines for Fortran statements

1. Fortran2003 deprecated features will not be accepted in new code (equivalence, common blocks, block data, etc.).
2. Use statements – should nearly always use the **only** construct to define what variables are being used from a given module. The only case when this is not to be done is when the list of variables is so long to make this unmanageable taking up more than 6 continue lines.
3. The "intent" attribute of subroutine input/output variables is used. Argument lists should ALWAYS use the INTENT attribute (try to have data be either IN or OUT and avoid INOUT if possible/reasonable).
4. Try to limit the use of the "use statement" to other modules to subroutines/functions and read-only parameters.
5. Data that is read-only, or output from a subroutine/function should be passed down using argument lists.
6. Modules explicitly declare public subroutines/functions at the module level.
7. By default all modules should have variables/methods/types declared to be private.
8. Print statements are not allowed, use write instead.
9. Stop statements will not be accepted in new code – use the endrun subroutine part of the CLM abortutils module.
10. Debugging statements should be removed prior to check in. If there is a need to leave them in, then they will be controlled by a logical or integer control variable.
11. Don't hard-code "magic numbers". For example, hard-coding a variable to "5" in order to make sure a certain path of an if branch is taken. Use the parameter statement at the top of the subroutine to set such variables.
12. Use symbolic logical operators (\<, \>, == etc.) rather than word (.lt., .gt. .eq.).

### Guidelines for working within the CLM codebase

1. We no longer accept new code with preprocessor macros (CPP tokens).
2. Namelist items should be defined in the XML files in `bld/namelist_files/`.
3. Type declarations make use of the shr_kinds_mod module.
4. All write statements are to iulog from clm_varctl.F90 or to another appropriate logical unit.
5. Write statements should be within "if ( masterproc )" if statement unless output by all tasks is really desired (such as for error statements).
6. Unit numbers are to be controlled using the CLM module fileutils.F90.
7. The module spmdMod should be used to get MPI information such as masterproc, iam etc.
8. The ncdio_pio module should be used to read/write/define NetCDF input/output data.
9. Modules clmvarcon or shr_const_mod should be used to define parameters in common with CLM and/or CESM.
10. Module clm_time_manager should be used to get time/date information. In particular length of day, month, or year information should be queried from clm_time_manager.
11. CLM parameters should be retrieved from the clm_varcon, clm_varpar and clm_varctl modules.
12. If writing is done in an initialization subroutine iulog should be flushed using shr_sys_flush at the end of the subroutine.
13. Evaluate new history fields being added. If they are not critical leave out of the code base; if they are only needed rarely, make them "inactive" by default.

### Guidelines for adding new datasets

1. Input datasets or raw datasets used to create input datasets (such as those that require the mksurfdata tool) need to be imported into by one of the CLM maintainers, who have the necessary permissions to do this and can also stage data on DIN_LOC_ROOT for NCAR internal machines.
2. Input datasets should have names that describe what they are, contain their resolution, and a have creation date at the end (i.e. \_cYYMMDD). Data should be in NetCDF format ending in a .nc extension.
3. All input datasets will be entered into the CLM XML database in the `bld/namelist_files/` directory.
4. CF-metadata conventions should be followed for the files.
5. The attributes of "units" and "long_name" should be set for all variables. Variables with missing data should have set the "missing_value" and "\_FillValue" attributes.
6. The attributes "Conventions", "source" and/or "history" should be set as global attributes on the file (references, institution, title and date should be considered as well).
7. Filenames should not be hardcoded into CLM code, they should be entered in by a namelist.
8. File sizes or number of time samples or records, should not be hardcoded into CLM code, but read in from the file.
