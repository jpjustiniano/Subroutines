SUBROUTINE diajuliano (day, month, year, dayj)   
!This program calculates the day of year corresponding to a specified date.
implicit none

real, intent(in):: day          !Day (dd)
real, intent(in) :: month        !Month (mm)
real, intent(in) :: year         !Year (yyyy)
real, intent(out) :: dayj 		!Day of year
integer :: i            			!Index,variable
integer :: leap_day     			!Extra day for leap year
 
! Check for leap year, and add extra day if necessary
IF ( mod(year,400.) == 0 ) THEN
    leap_day = 1    ! Years divisible by 400 are leap years
ELSE IF ( mod(year,100.) == 0 ) THEN
    leap_day = 0    ! Other centuries are not leap years
ELSE IF ( mod(year,4.) == 0 ) THEN
    leap_day = 1    ! Otherwise every 4th year 1S a leap year
ELSE
    leap_day = 0    ! Other years are not leap years
END IF


! Calculate day of year
dayj= day
DO i = 1, nint(month)-1
    ! Add days in months from January to last month
    SELECT CASE (i)
    CASE (1,3,5,7,8,10,12)
    dayj = dayj + 31
    CASE (4,6,9,11)
    dayj = dayj + 30
    CASE (2)
    dayj = dayj + 28 + leap_day
    END SELECT
END DO

END SUBROUTINE diajuliano
