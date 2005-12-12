subroutine initstats
  use precision
  use stats

  implicit none

  npos = 0
  ndif = 0
  difmx = 0.
  rdifmx = 0.
  xmx(1) = -9.d88
  xmx(2) = -9.d88
  xmn(1) = +9.d88
  xmn(2) = +9.d88
  rms = 0.
  xbar(1) = 0.
  xbar(2) = 0.
  rdbar = 0.
  rdlnbar = 0.
  diffmw = 0.
  return
end subroutine initstats
