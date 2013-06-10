module dynlandMod

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: dynlandMod
!
! !USES:
   use spmdMod
   use clmtype
   use decompMod   , only : get_proc_bounds
   use clm_varctl  , only : iulog
   use shr_kind_mod, only : r8 => shr_kind_r8
   use abortutils  , only : endrun
!
! !DESCRIPTION:
! Compute heat and water content to track conservation wrt dynamic land use
!
! !PUBLIC TYPES:
   implicit none
   private
   save
   public :: dynland_hwcontent
!
! !REVISION HISTORY:
!    2009-feb-20 B. Kauffman, created by
!
!EOP
!
! ! PRIVATE TYPES

!===============================================================================

contains
  
!===============================================================================
!BOP
!
! !ROUTINE: dynland_hwcontent
!
! !INTERFACE:

   subroutine dynland_hwcontent(begg,endg,gcell_liq,gcell_ice,gcell_heat)
 
! !DESCRIPTION:
!    Compute grid-level heat and water content
!
! !REVISION HISTORY:
!    2009-feb-20 B. Kauffman, created by
!
! !USES:

   use clm_varcon, only : istsoil,istice,istwet,istdlak,istslak,isturb,istice_mec
   use clm_varcon, only : istcrop
   use clm_varcon, only : icol_road_perv,icol_road_imperv,icol_roof
   use clm_varcon, only : icol_sunwall,icol_shadewall
   use clm_varcon, only : cpice,  cpliq
   use clm_varpar, only : nlevsno, nlevgrnd, nlevurb
   use clm_varpar, only : nlevlak
   use clm_varcon, only : denh2o

   implicit none

! !ARGUMENTS:

   integer , intent(in)  :: begg, endg              ! proc beg & end gridcell indices
   real(r8), intent(out) :: gcell_liq(begg:endg)
   real(r8), intent(out) :: gcell_ice  (begg:endg)
   real(r8), intent(out) :: gcell_heat (begg:endg)
 
! !LOCAL VARIABLES:
!EOP

   integer  :: li,lf         ! loop initial/final indicies
   integer  :: ci,cf         ! loop initial/final indicies
   integer  :: pi,pf         ! loop initial/final indicies

   integer  :: g,l,c,p,k     ! loop indicies (grid,lunit,column,pft,vertical level)

   real(r8) :: wtgcell       ! weight relative to grid cell
   real(r8) :: wtcol         ! weight relative to column
   real(r8) :: liq           ! sum of liquid water at column level
   real(r8) :: ice           ! sum of frozen water at column level
   real(r8) :: heat          ! sum of heat content at column level
   real(r8) :: cv            ! heat capacity [J/(m^2 K)]

   logical ,pointer :: pactive(:)        ! true=>do computations on this pft (see reweightMod for details)
    
   integer ,pointer :: ltype(:)          ! landunit type index
   integer ,pointer :: ctype(:)          ! column   type index
   integer ,pointer :: ptype(:)          ! pft      type index

   integer,  pointer :: nlev_improad(:)  ! number of impervious road layers
   real(r8), pointer :: cv_wall(:,:)     ! thermal conductivity of urban wall
   real(r8), pointer :: cv_roof(:,:)     ! thermal conductivity of urban roof
   real(r8), pointer :: cv_improad(:,:)  ! thermal conductivity of urban impervious road

   integer , pointer :: snl(:)           ! number of snow layers
   real(r8), pointer :: t_soisno(:,:)    ! soil temperature (Kelvin)
   real(r8), pointer :: h2osno(:)        ! snow water (mm H2O)
   real(r8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
   real(r8), pointer :: h2osoi_ice(:,:)  ! frozen water (kg/m2)
   real(r8), pointer :: watsat(:,:)      ! volumetric soil water at saturation (porosity)
   real(r8), pointer :: csol(:,:)        ! heat capacity, soil solids (J/m**3/Kelvin)
   real(r8), pointer :: dz(:,:)          ! layer depth (m)
   real(r8), pointer :: wa(:,:)          ! h2o in underground aquifer


!-------------------------------------------------------------------------------
! Note: this routine does not compute heat or water content of lakes.
!
!-------------------------------------------------------------------------------

   ! Set pointers into derived type


   pactive => pft%active

   ltype => lun%itype
   ctype => col%itype
   ptype =>pft%itype

   nlev_improad => lps%nlev_improad
   cv_wall      => lps%cv_wall
   cv_roof      => lps%cv_roof
   cv_improad   => lps%cv_improad

   snl          => cps%snl
   watsat       => cps%watsat
   csol         => cps%csol
   dz           => cps%dz
   t_soisno     => ces%t_soisno
   h2osoi_liq   => cws%h2osoi_liq
   h2osoi_ice   => cws%h2osoi_ice
   h2osno       => cws%h2osno

   ! Get relevant sizes

   do g = begg,endg ! loop over grid cells

      gcell_liq  (g) = 0.0_r8   ! sum for one grid cell
      gcell_ice  (g) = 0.0_r8   ! sum for one grid cell
      gcell_heat (g) = 0.0_r8   ! sum for one grid cell

      li = grc%luni(g)
      lf = grc%lunf(g)
      do l = li,lf   ! loop over land units  

         ci = lun%coli(l)
         cf = lun%colf(l)
         do c = ci,cf   ! loop over columns

            liq   = 0.0_r8 ! sum for one column
            ice   = 0.0_r8
            heat  = 0.0_r8

            !--- water & ice, above ground only ---
            if ( (ltype(l) == istsoil .or. ltype(l) == istcrop         )  &
            .or. (ltype(l) == istwet                                   )  &
            .or. (ltype(l) == istice                                   )  &
            .or. (ltype(l) == istice_mec                               )  &           
            .or. (ltype(l) == isturb .and. ctype(c) == icol_roof       )  &
            .or. (ltype(l) == isturb .and. ctype(c) == icol_road_imperv)  &
            .or. (ltype(l) == istdlak                                  )  &
            .or. (ltype(l) == isturb .and. ctype(c) == icol_road_perv  )) then

               if ( snl(c) < 0 ) then
                  do k = snl(c)+1,0 ! loop over snow layers
                     liq   = liq   + cws%h2osoi_liq(c,k)
                     ice   = ice   + cws%h2osoi_ice(c,k)
                  end do
               else                 ! no snow layers exist
                  ice = ice + cws%h2osno(c)
               end if
            end if

            !--- water & ice, below ground only ---
            if ( (ltype(l) == istsoil .or. ltype(l) == istcrop         )  &
            .or. (ltype(l) == istwet                                   )  &
            .or. (ltype(l) == istice                                   )  &
            .or. (ltype(l) == istdlak                                  )  &
            .or. (ltype(l) == istice_mec                               )  &           
            .or. (ltype(l) == isturb .and. ctype(c) == icol_road_perv  )) then
               do k = 1,nlevgrnd
                  liq   = liq   + cws%h2osoi_liq(c,k)
                  ice   = ice   + cws%h2osoi_ice(c,k)
               end do
            end if

            !--- water & ice, below ground, for lakes ---
            if ( ltype(l) == istdlak ) then
               do k = 1,nlevlak
                  liq   = liq   + (1 - cws%lake_icefrac(c,k))*cps%dz_lake(c,k)*denh2o
                  ice   = ice   + cws%lake_icefrac(c,k)*cps%dz_lake(c,k)*denh2o
                        ! lake layers do not change thickness when freezing, so denh2o should be used
                        ! (thermal properties are appropriately adjusted; see SLakeTemperatureMod)
               end do
            end if

            !--- water in aquifer ---
            if ( (ltype(l) == istsoil .or. ltype(l) == istcrop         )  &
            .or. (ltype(l) == istwet                                   )  &
            .or. (ltype(l) == istice                                   )  &
            .or. (ltype(l) == istice_mec                               )  &           
            .or. (ltype(l) == isturb .and. ctype(c) == icol_road_perv  )) then
               liq = liq + cws%wa(c)
            end if

            !--- water in canopy (at pft level) ---
            if (ltype(l) == istsoil .or. ltype(l) == istcrop) then   ! note: soil specified at LU level
               pi = col%pfti(c)
               pf = col%pftf(c)
               do p = pi,pf ! loop over pfts
                  if (pactive(p)) then
                     wtcol = pft%wtcol(p)
                     liq = liq + pws%h2ocan(p) * wtcol
                  end if
               end do
            end if

            if ( (ltype(l) /= istslak) ) then
            ! in new lake code, all lakes are "istdlak", but have variable depth

               !--- heat content, below ground only ---
               if (nlevurb > 0) then
                 do k = 1,nlevurb
                    if (ctype(c)==icol_sunwall .OR. ctype(c)==icol_shadewall) then
                        cv = cv_wall(l,k) * dz(c,k)
                        heat = heat + cv*t_soisno(c,k) / 1.e6_r8 
                    else if (ctype(c) == icol_roof) then
                        cv = cv_roof(l,k) * dz(c,k)
                        heat = heat + cv*t_soisno(c,k) / 1.e6_r8 
                    end if
                 end do
               end if
               do k = 1,nlevgrnd
                   if (ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall &
                       .and. ctype(c) /= icol_roof) then
                      if (ctype(c) == icol_road_imperv .and. k >= 1 .and. k <= nlev_improad(l)) then
                         cv = cv_improad(l,k) * dz(c,k)
                      else if (ltype(l) /= istwet .AND. ltype(l) /= istice .AND. ltype(l) /= istice_mec) then
                         cv = csol(c,k)*(1-watsat(c,k))*dz(c,k) + (h2osoi_ice(c,k)*cpice + h2osoi_liq(c,k)*cpliq)
                      else
                         cv = (h2osoi_ice(c,k)*cpice + h2osoi_liq(c,k)*cpliq)
                      endif
                      heat = heat + cv*t_soisno(c,k) / 1.e6_r8 
                   end if
               end do

               !--- heat content, below ground in lake water, for lakes ---
               do k = 1,nlevlak
                  if (ltype(l) == istdlak) then
                     cv = denh2o*cps%dz_lake(c,k)*( cws%lake_icefrac(c,k)*cpice + &
                          (1 - cws%lake_icefrac(c,k))*cpliq )
                     heat = heat + cv*ces%t_lake(c,k) / 1.e6_r8
                  end if
               end do

               !--- heat content, above ground only ---
               if ( snl(c) < 0 ) then
                  do k = snl(c)+1,0 ! loop over snow layers
                     cv = cpliq*h2osoi_liq(c,k) + cpice*h2osoi_ice(c,k)
                     heat = heat + cv*t_soisno(c,k) / 1.e6_r8
                  end do
               else if ( h2osno(c) > 0.0_r8 .and. ltype(l) /= istdlak) then
                       ! the heat capacity (not latent heat) of snow without snow layers
                       ! is currently ignored in SLakeTemperature, so it should be ignored here
                  k = 1
                  cv = cpice*h2osno(c)
                  heat = heat + cv*t_soisno(c,k) / 1.e6_r8
               end if

            end if

            !--- scale x/m^2 column-level values into x/m^2 gridcell-level values ---
            wtgcell = col%wtgcell(c)
            gcell_liq  (g) = gcell_liq  (g) + liq   * wtgcell
            gcell_ice  (g) = gcell_ice  (g) + ice   * wtgcell
            gcell_heat (g) = gcell_heat (g) + heat  * wtgcell

         end do ! column loop      
      end do ! landunit loop
   end do ! grid cell loop

   end subroutine dynland_hwcontent

!===============================================================================

end module dynlandMod
