module MOM_error_handler

  implicit none ; private

  public MOM_error, MOM_mesg, FATAL, WARNING

  integer, parameter :: FATAL = 0
  integer, parameter :: WARNING= 0

  contains

  subroutine MOM_error(level, message, all_print)
    integer,           intent(in) :: level
    character(len=*),  intent(in) :: message
    logical, optional, intent(in) :: all_print
  end subroutine MOM_error

  subroutine MOM_mesg(message, verb, all_print)
  character(len=*), intent(in)  :: message
  integer, optional, intent(in) :: verb
  logical, optional, intent(in) :: all_print
  ! This provides a convenient interface for writing an informative comment.
  integer :: verb_msg
  logical :: write_msg

end subroutine MOM_mesg

end module MOM_error_handler
