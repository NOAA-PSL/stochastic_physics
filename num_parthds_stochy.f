!>@brief The function 'num_parthds_stochy' sets the number of threads from the OMP_NUM_THREADS enviroment variable
!>@details This code is taken from the legacy spectral GFS
      function num_parthds_stochy()
      integer:: number_of_openMP_threads
      character(2) :: omp_threads
      integer :: stat
      call get_environment_variable("OMP_NUM_THREADS",omp_threads)
      read(omp_threads,*,iostat=stat)number_of_openMP_threads
      num_parthds_stochy = number_of_openMP_threads
      return
      end

