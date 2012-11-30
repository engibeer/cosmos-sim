
!	parameter ( Np_chunk = 1  000 000 )
!	integer ipc(Np_chunk)   ! parent cell index
	integer ipc(npmax)   ! parent cell index
        common / NEWDENS00 / ipc  
!        common / NEWDENS01 / dpc(mcell,Nproc)
!        common / NEWDENS02 / dpcR(mcell,Nproc)  
	integer iPL(npmax)  ! level of the cell to which particle belongs 
	common / PART21 / iPL   
