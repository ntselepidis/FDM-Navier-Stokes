program navier_stokes_simulation
  use m_sim
  use m_nvtx, only: nvtx_push, nvtx_pop
  implicit none
  type(simulation) :: sim
  real             :: time
  integer          :: step

  call sim%read_inputs()

  call sim%initialize()

  time = 0.0; step = 0
  do while ( time  < sim%total_time )
    ! evolve simulation
    call nvtx_push('timestep')
    call sim%evolve(time)
    call nvtx_pop()
    ! advance timestep
    time = time + sim%dt
    ! increment number of timesteps
    step = step + 1
    ! print info every 20 timesteps
    if (modulo(step, 20) == 0) print*, 'Time, step:', time, step
  end do

  print*, 'Time, step:', time, step

  call sim%export_results()

  call sim%destroy()

end program navier_stokes_simulation
