module MOM_hor_index

  implicit none ; private

!> Container for horizontal index ranges for data, computational and global domains
type, public :: hor_index_type
  integer :: isc !< The start i-index of cell centers within the computational domain
  integer :: iec !< The end i-index of cell centers within the computational domain
  integer :: jsc !< The start j-index of cell centers within the computational domain
  integer :: jec !< The end j-index of cell centers within the computational domain

  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain

  integer :: isg !< The start i-index of cell centers within the global domain
  integer :: ieg !< The end i-index of cell centers within the global domain
  integer :: jsg !< The start j-index of cell centers within the global domain
  integer :: jeg !< The end j-index of cell centers within the global domain

  integer :: IscB !< The start i-index of cell vertices within the computational domain
  integer :: IecB !< The end i-index of cell vertices within the computational domain
  integer :: JscB !< The start j-index of cell vertices within the computational domain
  integer :: JecB !< The end j-index of cell vertices within the computational domain

  integer :: IsdB !< The start i-index of cell vertices within the data domain
  integer :: IedB !< The end i-index of cell vertices within the data domain
  integer :: JsdB !< The start j-index of cell vertices within the data domain
  integer :: JedB !< The end j-index of cell vertices within the data domain

  integer :: IsgB !< The start i-index of cell vertices within the global domain
  integer :: IegB !< The end i-index of cell vertices within the global domain
  integer :: JsgB !< The start j-index of cell vertices within the global domain
  integer :: JegB !< The end j-index of cell vertices within the global domain

  integer :: idg_offset !< The offset between the corresponding global and local i-indices.
  integer :: jdg_offset !< The offset between the corresponding global and local j-indices.
  logical :: symmetric  !< True if symmetric memory is used.
end type hor_index_type

end module MOM_hor_index

