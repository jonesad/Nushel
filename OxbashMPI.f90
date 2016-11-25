! Fortran code for running oxbash with mpi
Program OxbashMPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine distributeLoad(nSize, nNames)
!Distribute the worklad of nNames calculations to the nSize processes 
Integer nSize, nNames, nIO, nI, nIdx, nJIdx
Real rTime
Integer Dimension (:,:), Allocatable :: aanIdxs
Integer Dimension (:), Allocatable :: anNames
Integer Dimension (:), Allocatable :: anCounts
Real Dimension (:) Allocatable :: arFtimes 
Real Dimension (:) Allocatable :: arFcopy 
Real Dimension (:) Allocatable :: arSums 
Logical lTimes
 
Allocate(anIdxs(nSize, Max(nNames-nSize, 1)))
Allocate(anCounts(nSize))
Allocate(arSums(nSize))
Allocate(anNames(nNames))
Allocate(arFtimes(nNames))
Allocate(arFcopy(nNames))

anIdxs(:,:) = 0
arSums(:) = 0 
anCounts(:) = 0

Inquire(FILE='times.dat', EXISTS=lTimes)
If (lTimes) Then 
  Open(UNIT=8, FILE='times.dat', ACTION='READ')
  nI=0
  Do
    nI+=1
    Read(8,*, IOSTAT=nIO) nIdx, rTime
    If(nIO>0) Then
      Write(*,*) 'Error: invalid input in times.dat'
    Else If (nIO<0) Then
      End Do
    Else
      anNames(nI)=nIdx
      arFtimes(nI)=rTime
      arFcopy(nI)=rTime
  Close(8)
  Do nI = 0, nNames
    nIdx = MINLOC(arSums)      
    nJIdx = MAXLOC(arFtimes)
    arSums(nIdx) = arSums(nIdx) + arFtimes(nJIdx)
    anCounts(nIdx) = anCounts(nIdx) + 1
    aanIdxs(nIdx, anCounts(nIdx)) = arFcopy(nJIdx)
    
  
