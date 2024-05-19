      program intake_MD_w_inlet

      implicit none
      integer i,j,k
      integer npart,nspart,nlpart,maxnpart
      parameter (maxnpart=4096)
      double precision x(maxnpart),y(maxnpart),vx(maxnpart),vy(maxnpart)
      double precision mass(maxnpart)
      double precision radius(maxnpart)
      integer eject_label(maxnpart) !0: not ejected; 1: ejected
      double precision boxx,boxy,phi
      character*60 str1,str2,str3,str4
      double precision totV
      !double precision old_totV,rel_totV,min_rel_totV
      double precision stop_min_totV
      double precision stop_min_tot_ke,tot_ke
      double precision pi
      parameter(pi=3.141592653589793d0)
      double precision Rratio,radius1
      double precision dt,dt_sigma1
      integer skip
      double precision rcut,rlist,ranrlist
      integer n_i,maxn_i ! n_i = num of inclusion
      parameter(maxn_i=10)
      double precision i_x(maxn_i),i_y(maxn_i)
      double precision i_vx(maxn_i),i_vy(maxn_i)
      double precision i_radius(maxn_i),i_mass(maxn_i)
      integer MDstep_idx,last_MDstep_idx
      integer movie_switch
      double precision gamman,gamman_sigma1
      double precision gammat,gammat_sigma1
      common /damping/ gamman,gamman_sigma1,gammat,gammat_sigma1
      double precision mu
      common /self_propel/ mu
      integer dir_switch !0: right; 1: up; 2: down
      common /motion_mode/ dir_switch
      integer num_bump,bump_idx,num_intake_part
      double precision intake_range,ejection_pos
      double precision ran2,gasdev
      integer seed
      double precision alpha_x,alpha_y
      double precision roving_dist_x,roving_dist_y

      write(*,*) 'Enter name of file containing relaxed IC'
      read(*,*) str1
      write(*,*) 'Enter name of file saving part ejection record'
      read(*,*) str2
      write(*,*) 'Save moive? (0: No; 1: IC+all; 2: IC+ejections only)'
      read(*,*) movie_switch
      if(movie_switch .ne. 0) then
        write(*,*) 'Enter name of file saving configurations'
        read(*,*) str3
      endif
      write(*,*) 'Enter name of file saving summary data'
      read(*,*) str4
      write(*,*) 'Enter large/small particle size = Rratio >= 1.0'
      read(*,*) Rratio
      write(*,*) 'Enter value of time step of particle size 1'
      read(*,*) dt_sigma1
      write(*,*) 'Enter gamman_sigma1'
      read(*,*) gamman_sigma1
      write(*,*) 'Enter gammat_sigma1'
      read(*,*) gammat_sigma1
      write(*,*) 'Enter mu of self propulsion'
      read(*,*) mu
      write(*,*) 'Enter value of rlist/rcut'
      read(*,*) ranrlist
      !write(*,*) 'Enter min_rel_totV'
      !write(*,*) 'If rel_totV < min_rel_totV, relaxation stops.'
      !read(*,*) min_rel_totV
      write(*,*) 'Enter value of stop_min_totV for ending relaxation'
      write(*,*) '(condition 1/2)'
      read(*,*) stop_min_totV
      write(*,*) 'Enter value of stop_min_tot_ke for ending relaxation'
      write(*,*) '(condition 2/2)'
      read(*,*) stop_min_tot_ke
      write(*,*) 'Enter skip for saving data'
      read(*,*) skip
      write(*,*) 'Enter moving direction (0: right; 1: up; 2: down)'
      read(*,*) dir_switch
      write(*,*) 'Enter num_bump for moving the rover'
      read(*,*) num_bump
      write(*,*) 'Enter intake_range in unit of radius1'
      read(*,*) intake_range
      write(*,*) 'Enter ejection_pos in unit of radius1'
      read(*,*) ejection_pos
      write(*,*) 'Enter alpha_x (>= 0) of Gaussian roving dist'
      write(*,*) 'in unit of 1d-2 * i_radius(1) in x direction'
      read(*,*) alpha_x
      write(*,*) 'Enter alpha_y (>= 0) of Gaussian roving dist'
      write(*,*) 'in unit of 1d-2 * i_radius(1) in y direction'
      read(*,*) alpha_y

      open(unit=1,file=str1,status='old')
      open(unit=2,file=str2,status='unknown')
      if(movie_switch .ne. 0) then
        open(unit=3,file=str3,status='unknown')
      endif
      open(unit=4,file=str4,status='unknown')
 40   format((I5,3X),(I3,3X),4(e21.14,3X),(I12,3X))

      read(1,*)
      read(1,*) nspart,nlpart,n_i,phi
      npart = nspart + nlpart

      do i=1,npart
        read(1,*) x(i), y(i), vx(i), vy(i), radius(i)
      enddo

      do i=1,n_i
        read(1,*) i_x(i),i_y(i),i_vx(i),i_vy(i),i_radius(i)
      enddo

      boxx = 1.d0
      boxy = 1.d0

      radius1 = radius(1)
      write(*,*) 'radius1 =',radius1
      rcut = 2.d0 * radius1 * Rratio
      write(*,*) 'rcut =',rcut
      rlist = ranrlist * rcut

      !----------------------   setting parameters   -------------------
      ! the following lines translate length scale from particle size 1
      ! to box size 1. Choose smallest particle diameter as the ruler.
      dt = dt_sigma1 * 2.d0 * radius1
      gamman = gamman_sigma1 / (2.d0 * radius1)
      gammat = gammat_sigma1 / (2.d0 * radius1)
      !-----------------------------------------------------------------

      write(2,*) '# bump_idx, part_idx (part election record)'
      write(4,*) '# bump_idx, #_intake_part, i_x, i_y, i_vx, i_vy, seed'

      do i=1,npart
        mass(i) = (radius(i) / radius1) * (radius(i) / radius1)
      enddo

      do i=1,n_i
        i_mass(i) = (i_radius(i) / radius1) * (i_radius(i) / radius1)
      enddo



      do i=1,npart
        eject_label(i) = 0 !0: not ejected
      enddo

      ! save IC for movie
 30   format(5(e21.14,3X),I1)
 31   format(5(e21.14,3X))

      if(movie_switch .ne. 0) then
        write(3,*) npart + n_i, 0 ! bump_idx = 0 = IC
        write(3,*) nspart,nlpart,n_i,phi
        do i=1,npart
          write(3,30) x(i),y(i),vx(i),vy(i),radius(i),eject_label(i)
        enddo
        do i=1,n_i
          write(3,31) i_x(i),i_y(i),i_vx(i),i_vy(i),i_radius(i)
        enddo
        write(3,*) ''
        backspace(3)
        read(3,*)
      endif

      ! initialization of random number seed
      seed = -1

      ! initialization
      last_MDstep_idx = 0

      do bump_idx=0,num_bump

        ! when bump_idx = 0, relax IC only.

        if(bump_idx .gt. 0) then
          !------------------------- rover bumping -------------------------
          num_intake_part = 0

          if(dir_switch .eq. 0) then !0: roving right
            do i=1,npart
              if(    x(i) .ge. i_x(1) + i_radius(1)-intake_range*radius1
     &         .and. x(i) .le. i_x(1) + i_radius(1)
     &         .and. y(i) .ge. i_y(1) - intake_range/2.d0*radius1 
     &         .and. y(i) .le. i_y(1) + intake_range/2.d0*radius1) then
                x(i) = x(i) - i_radius(1)
     &                      - ejection_pos * radius1
                num_intake_part = num_intake_part + 1
                eject_label(i) = 1 !1: ejected

                write(2,*) bump_idx,i
                backspace(2)
                read(2,*)
              endif
            enddo
          elseif(dir_switch .eq. 1) then !1: roving up
            do i=1,npart
              if(    y(i) .ge. i_y(1) + i_radius(1)-intake_range*radius1
     &         .and. y(i) .le. i_y(1) + i_radius(1)
     &         .and. x(i) .ge. i_x(1) - intake_range/2.d0*radius1
     &         .and. x(i) .le. i_x(1) + intake_range/2.d0*radius1) then
                y(i) = y(i) - i_radius(1)
     &                      - ejection_pos * radius1
                num_intake_part = num_intake_part + 1
                eject_label(i) = 1 !1: ejected

                write(2,*) bump_idx,i
                backspace(2)
                read(2,*)
              endif
            enddo
          elseif(dir_switch .eq. 2) then !2: roving down
            do i=1,npart
              if(    y(i) .le. i_y(1) - i_radius(1)+intake_range*radius1
     &         .and. y(i) .ge. i_y(1) - i_radius(1)
     &         .and. x(i) .ge. i_x(1) - intake_range/2.d0*radius1 
     &         .and. x(i) .le. i_x(1) + intake_range/2.d0*radius1) then
                y(i) = y(i) + i_radius(1)
     &                      + ejection_pos * radius1
                num_intake_part = num_intake_part + 1
                eject_label(i) = 1 !1: ejected

                write(2,*) bump_idx,i
                backspace(2)
                read(2,*)
              endif
            enddo
          endif

          write(*,*) '#',bump_idx,'th bump:'
          write(*,*) '# num_intake_part =',num_intake_part

          if(num_intake_part .eq. 0) then
            ! move rover in the dir_switch direction !0: right; 1: up; 2: down
            if(dir_switch .eq. 0) then !0: roving right
              roving_dist_x = dabs(gasdev(seed)) * alpha_x
     &                      * 1.d-2 * i_radius(1)
              roving_dist_y = gasdev(seed) * alpha_y
     &                      * 1.d-2 * i_radius(1)
            elseif(dir_switch .eq. 1) then !1: roving up
              roving_dist_x = gasdev(seed) * alpha_x
     &                      * 1.d-2 * i_radius(1)
              roving_dist_y = dabs(gasdev(seed)) * alpha_y
     &                      * 1.d-2 * i_radius(1)
            elseif(dir_switch .eq. 2) then !2: roving down
              roving_dist_x = gasdev(seed) * alpha_x
     &                      * 1.d-2 * i_radius(1)
              roving_dist_y = - dabs(gasdev(seed)) * alpha_y
     &                      * 1.d-2 * i_radius(1)
            endif

            i_x(1) = i_x(1) + roving_dist_x
            i_y(1) = i_y(1) + roving_dist_y
          endif

          if(movie_switch .eq. 1) then
            ! save config after each bump
            write(3,*) npart + n_i, bump_idx
            write(3,*) nspart,nlpart,n_i,phi
            do i=1,npart
              write(3,30) x(i),y(i),vx(i),vy(i),radius(i),eject_label(i)
            enddo
            do i=1,n_i
              write(3,31) i_x(i),i_y(i),i_vx(i),i_vy(i),i_radius(i)
            enddo
            write(3,*) ''
            backspace(3)
            read(3,*)
          elseif(movie_switch .eq. 2) then
            ! save config after each bump w/ nonzero num_intake_part
            if(num_intake_part .ne. 0) then
              write(3,*) npart + n_i, bump_idx
              write(3,*) nspart,nlpart,n_i,phi
              do i=1,npart
                write(3,30) x(i),y(i),vx(i),vy(i),radius(i),
     &                      eject_label(i)
              enddo
              do i=1,n_i
                write(3,31) i_x(i),i_y(i),i_vx(i),i_vy(i),i_radius(i)
              enddo
              write(3,*) ''
              backspace(3)
              read(3,*)
            endif
          endif
          !------------------------- rover bumping -------------------------
        endif !if(bump_idx .gt. 0) then

        !------------------------ MD relaxation ------------------------
        !initiallize
        MDstep_idx = 0
        !old_totV = 1.d0  !a nonzero number
        !rel_totV = 1.d0  !a number > min_rel_totV

        !do while(rel_totV .ge. min_rel_totV)
 1      continue
          call verlet(x,y,i_x,i_y,vx,vy,i_vx,i_vy,npart,n_i,boxx,boxy,
     &                radius,i_radius,dt,rcut,rlist,mass,i_mass,totV,
     &                intake_range,radius1)

          MDstep_idx = MDstep_idx + 1

          call KE_cal(vx,vy,mass,i_vx,i_vy,i_mass,tot_ke,npart,n_i)

          if(mod(MDstep_idx,skip) .eq. 0) then

            !rel_totV = dabs(dabs(totV) - old_totV) / old_totV
            !old_totV = dabs(totV)

            write(*,*) 'step',MDstep_idx,'t* =',MDstep_idx * dt_sigma1
            write(*,*) 'totV =',totV
            write(*,*) 'tot_ke =',tot_ke
            !write(*,*) 'rel_totV =',rel_totV
            write(*,*) ' '

          endif !if(mod(MDstep_idx,skip) .eq. 0)

          if(totV   .lt. stop_min_totV .and. 
     &       tot_ke .lt. stop_min_tot_ke) goto 10
        goto 1
 10     continue
        !enddo !do while(rel_totV .ge. min_rel_totV)

        write(*,*) 'stop at step', MDstep_idx
        last_MDstep_idx = last_MDstep_idx + MDstep_idx

        write(*,*) 'relaxation stopped due to:'
        write(*,*) 'totV =',totV,'<',stop_min_totV,'and'
        write(*,*) 'tot_ke =',tot_ke,'<',stop_min_tot_ke
        write(*,*) ''

        if(movie_switch .eq. 1) then
          ! save the IC and all configs
          write(3,*) npart + n_i, bump_idx
          write(3,*) nspart,nlpart,n_i,phi
          do i=1,npart
            write(3,30) x(i),y(i),vx(i),vy(i),radius(i),eject_label(i)
          enddo
          do i=1,n_i
            write(3,31) i_x(i),i_y(i),i_vx(i),i_vy(i),i_radius(i)
          enddo
          write(3,*) ''
          backspace(3)
          read(3,*)
        elseif(movie_switch .eq. 2) then
          ! save only the IC and other configs involving ejection
          if(bump_idx .eq. 0 .or. num_intake_part .ne. 0) then
            write(3,*) npart + n_i, bump_idx
            write(3,*) nspart,nlpart,n_i,phi
            do i=1,npart
              write(3,30) x(i),y(i),vx(i),vy(i),radius(i),eject_label(i)
            enddo
            do i=1,n_i
              write(3,31) i_x(i),i_y(i),i_vx(i),i_vy(i),i_radius(i)
            enddo
            write(3,*) ''
            backspace(3)
            read(3,*)
          endif
        endif

        if(bump_idx .eq. 0) then
          write(4,40) bump_idx, 0,
     &                i_x(1), i_y(1), i_vx(1), i_vy(1), seed
        else
          write(4,40) bump_idx, num_intake_part,
     &                i_x(1), i_y(1), i_vx(1), i_vy(1), seed
        endif
        !------------------------ MD relaxation ------------------------

        !-- a protection because there is no cruiser-box interaction ---
        if(i_x(1) + i_radius(1) .gt.  0.5d0 * boxx .or.
     &     i_x(1) - i_radius(1) .lt. -0.5d0 * boxx .or.
     &     i_y(1) + i_radius(1) .gt.  0.5d0 * boxy .or.
     &     i_y(1) - i_radius(1) .lt. -0.5d0 * boxy) then
          write(*,*) 'The cruiser is outside the simulation box, stop!'
          stop
        endif
        !-- a protection because there is no cruiser-box interaction ---

      enddo !do bump_idx=0,num_bump

      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine verlet(x,y,i_x,i_y,vx,vy,i_vx,i_vy,npart,n_i,boxx,boxy,
     &                  radius,i_radius,dt,rcut,rlist,mass,i_mass,totV,
     &                  intake_range,radius1)

      ! See Allen & Tildesley p.81, Eq.(3.19)
      implicit none
      integer i,j,k
      integer npart,maxnpart
      parameter (maxnpart=4096)
      double precision x(maxnpart),y(maxnpart),vx(maxnpart),vy(maxnpart)
      double precision radius(maxnpart)
      double precision boxx,boxy
      double precision dt
      double precision rcut,rlist
      double precision fx(maxnpart),fy(maxnpart),V
      double precision dampfnx(maxnpart),dampfny(maxnpart)
      double precision dampftx(maxnpart),dampfty(maxnpart)
      double precision wallfx(maxnpart),wallfy(maxnpart),wallV
      double precision inletwallfx(maxnpart),inletwallfy(maxnpart)
      double precision inletwallV
      double precision mass(maxnpart),G_fy(maxnpart)
      double precision self_propel_fx(maxnpart),self_propel_fy(maxnpart)
      double precision old_ax(maxnpart),old_ay(maxnpart)
      double precision new_ax(maxnpart),new_ay(maxnpart)
      double precision inclusionfx(maxnpart),inclusionfy(maxnpart)
      integer n_i,maxn_i ! n_i = num of inclusion
      parameter(maxn_i=10)
      double precision i_x(maxn_i),i_y(maxn_i)
      double precision i_fx(maxn_i),i_fy(maxn_i)
      double precision i_inletwallfx(maxn_i),i_inletwallfy(maxn_i)
      double precision i_vx(maxn_i),i_vy(maxn_i)
      double precision i_radius(maxn_i),i_mass(maxn_i),i_G_fy(maxn_i)
      double precision i_self_propel_fx(maxn_i),i_self_propel_fy(maxn_i)
      double precision old_i_ax(maxn_i),old_i_ay(maxn_i)
      double precision new_i_ax(maxn_i),new_i_ay(maxn_i)
      double precision G_V,totV
      double precision intake_range,radius1

      do i=1,npart
        old_ax(i) = 0.d0
        old_ay(i) = 0.d0
        new_ax(i) = 0.d0
        new_ay(i) = 0.d0
      enddo

      do i=1,n_i
        old_i_ax(i) = 0.d0
        old_i_ay(i) = 0.d0
        new_i_ax(i) = 0.d0
        new_i_ay(i) = 0.d0
      enddo

      call Vlistforce(rcut,rlist,x,y,fx,fy,V,npart,boxx,boxy,radius)
      call Vlist_dampforce(rcut,rlist,x,y,vx,vy,dampfnx,dampfny,
     &                     dampftx,dampfty,npart,boxx,boxy,radius)
      call wallforce(x,y,wallfx,wallfy,wallV,npart,boxx,boxy,radius)
      call inclusion_force(x,y,i_x,i_y,inclusionfx,inclusionfy,
     &                     i_fx,i_fy,npart,n_i,
     &                     boxx,boxy,radius,i_radius,
     &                     intake_range,radius1)
      call inlet_wallforce(x,y,inletwallfx,inletwallfy,inletwallV,
     &                     i_x,i_y,i_inletwallfx,i_inletwallfy,
     &                     radius,npart,i_radius,n_i,
     &                     intake_range,radius1)
      call self_propel_force(vx,vy,self_propel_fx,self_propel_fy,npart)
      call inclusion_self_propel_force(i_vx,i_vy,
     &                                 i_self_propel_fx,
     &                                 i_self_propel_fy,n_i)

      do i=1,npart
        old_ax(i) = ( fx(i) + dampfnx(i) + dampftx(i)
     &            + wallfx(i) + inclusionfx(i)
     &            + inletwallfx(i) + self_propel_fx(i) )
     &            / mass(i)
        old_ay(i) = ( fy(i) + dampfny(i) + dampfty(i)
     &            + wallfy(i) + inclusionfy(i)
     &            + inletwallfy(i) + self_propel_fy(i) )
     &            / mass(i)
      enddo

      do i=1,npart
        x(i) = x(i) + vx(i) * dt + (old_ax(i) * dt * dt) / 2.d0
        y(i) = y(i) + vy(i) * dt + (old_ay(i) * dt * dt) / 2.d0
      enddo

      do i=1,n_i
        old_i_ax(i) = ( i_fx(i)
     &              + i_inletwallfx(i) + i_self_propel_fx(i) )
     &              / i_mass(i)
        old_i_ay(i) = ( i_fy(i) 
     &              + i_inletwallfy(i) + i_self_propel_fy(i) )
     &              / i_mass(i)
      enddo

      do i=1,n_i
        i_x(i) = i_x(i) + i_vx(i) * dt + (old_i_ax(i) * dt * dt) / 2.d0
        i_y(i) = i_y(i) + i_vy(i) * dt + (old_i_ay(i) * dt * dt) / 2.d0
      enddo

      call Vlistforce(rcut,rlist,x,y,fx,fy,V,npart,boxx,boxy,radius)
      call Vlist_dampforce(rcut,rlist,x,y,vx,vy,dampfnx,dampfny,
     &                     dampftx,dampfty,npart,boxx,boxy,radius)
      call wallforce(x,y,wallfx,wallfy,wallV,npart,boxx,boxy,radius)
      call inclusion_force(x,y,i_x,i_y,inclusionfx,inclusionfy,
     &                     i_fx,i_fy,npart,n_i,
     &                     boxx,boxy,radius,i_radius,
     &                     intake_range,radius1)
      call inlet_wallforce(x,y,inletwallfx,inletwallfy,inletwallV,
     &                     i_x,i_y,i_inletwallfx,i_inletwallfy,
     &                     radius,npart,i_radius,n_i,
     &                     intake_range,radius1)
      call self_propel_force(vx,vy,self_propel_fx,self_propel_fy,npart)
      call inclusion_self_propel_force(i_vx,i_vy,
     &                                 i_self_propel_fx,
     &                                 i_self_propel_fy,n_i)

      do i=1,npart
        new_ax(i) = ( fx(i) + dampfnx(i) + dampftx(i)
     &            + wallfx(i) + inclusionfx(i)
     &            + inletwallfx(i) + self_propel_fx(i) )
     &            / mass(i)
        new_ay(i) = ( fy(i) + dampfny(i) + dampfty(i)
     &            + wallfy(i) + inclusionfy(i)
     &            + inletwallfy(i) + self_propel_fy(i) )
     &            / mass(i)
      enddo

      do i=1,npart
        vx(i) = vx(i) + (new_ax(i) + old_ax(i)) / 2.d0 * dt
        vy(i) = vy(i) + (new_ay(i) + old_ay(i)) / 2.d0 * dt
      enddo

      do i=1,n_i
        new_i_ax(i) = ( i_fx(i)
     &              + i_inletwallfx(i) + i_self_propel_fx(i) )
     &              / i_mass(i)
        new_i_ay(i) = ( i_fy(i)
     &              + i_inletwallfy(i) + i_self_propel_fy(i) )
     &              / i_mass(i)
      enddo

      do i=1,n_i
        i_vx(i) = i_vx(i) + (new_i_ax(i) + old_i_ax(i)) / 2.d0 * dt
        i_vy(i) = i_vy(i) + (new_i_ay(i) + old_i_ay(i)) / 2.d0 * dt
      enddo

      totV = V + wallV + inletwallV

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check(rcut,rlist,x,y,x0,y0,npart,update)

      implicit none
      integer i,j,k
      integer npart,maxnpart
      parameter(maxnpart=4096)
      double precision x(maxnpart),y(maxnpart)
      double precision x0(maxnpart),y0(maxnpart)
      double precision rcut,rlist
      integer update
      double precision dispmax
      save

      dispmax = 0.d0

      do i=1,npart

        dispmax = max(dabs(x(i) - x0(i)),dispmax)
        dispmax = max(dabs(y(i) - y0(i)),dispmax)
        
      enddo

      dispmax = 2.d0 * 1.414213562373095d0 * dispmax

      if (dispmax .gt. (rlist-rcut)) then
        update = 1 ! list need be updated
      else
        update = 0 ! list need not be updated
      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Vlistforce(rcut,rlist,x,y,fx,fy,
     &                      V,npart,boxx,boxy,radius)

      implicit none
      integer i,j,k
      integer npart,maxnpart,maxlist
      parameter(maxnpart=4096)
      parameter(maxlist=50000) ! increase it with increasing npart
      double precision x(maxnpart),y(maxnpart),fx(maxnpart),fy(maxnpart)
      double precision V
      double precision boxx,boxy
      double precision radius(maxnpart)
      double precision rx,ry,r,r2,sigma,sigma2,ff
      double precision rcut,rlist,rlist2
      double precision x0(maxnpart),y0(maxnpart)
      integer update,nlist,point(maxlist),list(maxlist),jbegin,jend,jn

      rlist2 = rlist * rlist

      V = 0.d0
      
      do i=1,npart
        fx(i) = 0.d0
        fy(i) = 0.d0
      enddo

      call check(rcut,rlist,x,y,x0,y0,npart,update)

      if (update .eq. 1) then ! list need be updated

        do i=1,npart
          x0(i) = x(i)
          y0(i) = y(i)
        enddo

        nlist = 0
       
        do i=1,npart-1
          point(i) = nlist + 1
          do j=i+1,npart

            rx = x(i) - x(j)
            ry = y(i) - y(j)
            !rx = rx - boxx*dnint(rx/boxx)
            !ry = ry - boxy*dnint(ry/boxy)
            r2 = (rx * rx + ry * ry)
          
            sigma = radius(i) + radius(j)
            sigma2 = sigma * sigma

            if (r2 .lt. rlist2) then

              nlist = nlist + 1
              list(nlist) = j

              if(r2 .lt. sigma2) then

                r = dsqrt(r2)

                ff = (1.d0 - r/sigma)

                fx(i) = fx(i) + (rx/r)*ff/sigma
                fy(i) = fy(i) + (ry/r)*ff/sigma
                fx(j) = fx(j) - (rx/r)*ff/sigma
                fy(j) = fy(j) - (ry/r)*ff/sigma

                V = V + 0.5d0 * ff * ff

              endif

            endif

          enddo
        enddo
        point(npart) = nlist + 1

      else ! list need not be updated

        do i=1,npart-1
          jbegin = point(i)
          jend = point(i+1) - 1
          if(jbegin .le. jend) then
            do jn=jbegin,jend
              j = list(jn)

              rx = x(i) - x(j)
              ry = y(i) - y(j)
              !rx = rx - boxx*dnint(rx/boxx)
              !ry = ry - boxy*dnint(ry/boxy)
              r2 = (rx*rx + ry*ry)

              sigma = radius(i) + radius(j)
              sigma2 = sigma*sigma

              if(r2 .lt. sigma2) then

                r = dsqrt(r2)

                ff = (1.d0 - r/sigma)

                fx(i) = fx(i) + (rx/r)*ff/sigma
                fy(i) = fy(i) + (ry/r)*ff/sigma
                fx(j) = fx(j) - (rx/r)*ff/sigma
                fy(j) = fy(j) - (ry/r)*ff/sigma

                V = V + 0.5d0*ff*ff

              endif
            enddo
          endif
        enddo

      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Vlist_dampforce(rcut,rlist,x,y,vx,vy,dampfnx,dampfny,
     &                           dampftx,dampfty,npart,boxx,boxy,radius)

      implicit none
      integer i,j,k
      integer npart,maxnpart,maxlist
      parameter(maxnpart=4096)
      parameter(maxlist=50000) ! increase it with increasing npart
      double precision x(maxnpart),y(maxnpart)
      double precision dampfnx(maxnpart),dampfny(maxnpart)
      double precision dampftx(maxnpart),dampfty(maxnpart)
      double precision vx(maxnpart),vy(maxnpart)
      double precision boxx,boxy,radius(maxnpart)
      double precision rx,ry,r,r2,sigma,sigma2
      double precision vxij,vyij
      double precision rcut,rlist,rlist2
      double precision x0(maxnpart),y0(maxnpart)
      integer update,nlist,point(maxlist),list(maxlist),jbegin,jend,jn
      double precision vn
      double precision gamman,gamman_sigma1
      double precision gammat,gammat_sigma1
      common /damping/ gamman,gamman_sigma1,gammat,gammat_sigma1
      save

      !----------------------   setting parameters   -------------------
      ! the following lines translate length scale from particle size 1
      ! to box size 1. Choose smaller particle diameter as the ruler.
      ! gamman = gamman_sigma1 / (2.d0 * radius(1))
      ! gammat = gammat_sigma1 / (2.d0 * radius(1))
      !-----------------------------------------------------------------

      rlist2 = rlist*rlist

      do i=1,npart
        dampfnx(i) = 0.d0
        dampfny(i) = 0.d0

        dampftx(i) = 0.d0
        dampfty(i) = 0.d0
      enddo

      call check(rcut,rlist,x,y,x0,y0,npart,update)

      if (update .eq. 1) then ! list need be updated

        do i=1,npart
          x0(i) = x(i)
          y0(i) = y(i)
        enddo

        nlist = 0

        do i=1,npart-1
          point(i) = nlist + 1
          do j=i+1,npart

            rx = x(i) - x(j)
            ry = y(i) - y(j)
            !rx = rx - boxx*dnint(rx/boxx)
            !ry = ry - boxy*dnint(ry/boxy)
            r2 = (rx*rx + ry*ry)

            vxij = vx(i) - vx(j)
            vyij = vy(i) - vy(j)

            sigma = radius(i) + radius(j)
            sigma2 = sigma*sigma

            if (r2 .lt. rlist2) then

              nlist = nlist + 1
              list(nlist) = j

              if(r2.lt.sigma2) then !only contact produces damping force

                r = dsqrt(r2) ! square root is an expensive operation

                vn = (vxij*rx + vyij*ry)/r

                !normal damping force is minus sign
                dampfnx(i) = dampfnx(i) - vn*(rx/r)*gamman
                dampfny(i) = dampfny(i) - vn*(ry/r)*gamman
                dampfnx(j) = dampfnx(j) + vn*(rx/r)*gamman
                dampfny(j) = dampfny(j) + vn*(ry/r)*gamman

                !tangential damping force is minus sign
                dampftx(i) = dampftx(i) - (vxij - vn*(rx/r))*gammat
                dampfty(i) = dampfty(i) - (vyij - vn*(ry/r))*gammat
                dampftx(j) = dampftx(j) + (vxij - vn*(rx/r))*gammat
                dampfty(j) = dampfty(j) + (vyij - vn*(ry/r))*gammat

              endif

            endif

          enddo
        enddo
        point(npart) = nlist + 1

      else ! list need not be updated

        do i=1,npart-1
          jbegin = point(i)
          jend = point(i+1) - 1
          if(jbegin .le. jend) then
            do jn=jbegin,jend
              j = list(jn)

              rx = x(i) - x(j)
              ry = y(i) - y(j)
              !rx = rx - boxx*dnint(rx/boxx)
              !ry = ry - boxy*dnint(ry/boxy)
              r2 = (rx*rx + ry*ry)

              vxij = vx(i) - vx(j)
              vyij = vy(i) - vy(j)

              sigma = radius(i) + radius(j)
              sigma2 = sigma*sigma

              if(r2.lt.sigma2) then !only contact produces damping force

                r = dsqrt(r2) ! square root is an expensive operation

                vn = (vxij*rx + vyij*ry)/r

                !damping force is minus sign
                dampfnx(i) = dampfnx(i) - vn*(rx/r)*gamman
                dampfny(i) = dampfny(i) - vn*(ry/r)*gamman
                dampfnx(j) = dampfnx(j) + vn*(rx/r)*gamman
                dampfny(j) = dampfny(j) + vn*(ry/r)*gamman

                !tangential damping force is minus sign
                dampftx(i) = dampftx(i) - (vxij - vn*(rx/r))*gammat
                dampfty(i) = dampfty(i) - (vyij - vn*(ry/r))*gammat
                dampftx(j) = dampftx(j) + (vxij - vn*(rx/r))*gammat
                dampfty(j) = dampfty(j) + (vyij - vn*(ry/r))*gammat

              endif
            enddo
          endif
       enddo

      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wallforce(x,y,wallfx,wallfy,wallV,
     &                     npart,boxx,boxy,radius)

      implicit none
      integer i,j,k
      integer npart,maxnpart
      parameter(maxnpart=4096)
      double precision x(maxnpart),y(maxnpart)
      double precision wallfx(maxnpart),wallfy(maxnpart)
      double precision imagex,imagey
      double precision boxx,boxy
      double precision wallV
      double precision radius(maxnpart)
      double precision rx,ry,r,r2,sigma,sigma2,ff

      wallV = 0.d0

      do i=1,npart
        wallfx(i) = 0.d0
        wallfy(i) = 0.d0
      enddo

      do i=1,npart

        ! hit left x-wall
        if( x(i) .lt. (-boxx/2.d0 + radius(i)) ) then
          imagex = -boxx - x(i)
          imagey = y(i)

          !__________________________________
          rx = x(i) - imagex
          ry = y(i) - imagey

          r2 = (rx**2 + ry**2)

          ! assume the particle and its image have the same size
          sigma = radius(i) + radius(i)

          r = dsqrt(r2)

          ff = (1.d0 - r/sigma)

          wallfx(i) = wallfx(i) + (rx/r)*ff/sigma
          wallfy(i) = wallfy(i) + (ry/r)*ff/sigma

          wallV = wallV + 0.5d0*ff**2
          !__________________________________
        endif

        ! hit right x-wall
        if( x(i) .gt. (boxx/2.d0 - radius(i)) ) then
          imagex = boxx - x(i)
          imagey = y(i)

          !__________________________________
          rx = x(i) - imagex
          ry = y(i) - imagey

          r2 = (rx**2 + ry**2)

          ! assume the particle and its image have the same size
          sigma = radius(i) + radius(i)

          r = dsqrt(r2)

          ff = (1.d0 - r/sigma)

          wallfx(i) = wallfx(i) + (rx/r)*ff/sigma
          wallfy(i) = wallfy(i) + (ry/r)*ff/sigma

          wallV = wallV + 0.5d0*ff**2
          !__________________________________
        endif

        ! hit lower y-wall
        if( y(i) .lt. (-boxy/2.d0 + radius(i)) ) then
          imagex = x(i)
          imagey = -boxy - y(i)

          !__________________________________
          rx = x(i) - imagex
          ry = y(i) - imagey

          r2 = (rx**2 + ry**2)

          ! assume the particle and its image have the same size
          sigma = radius(i) + radius(i)

          r = dsqrt(r2)

          ff = (1.d0 - r/sigma)

          wallfx(i) = wallfx(i) + (rx/r)*ff/sigma
          wallfy(i) = wallfy(i) + (ry/r)*ff/sigma

          wallV = wallV + 0.5d0*ff**2
          !__________________________________
        endif

        ! hit upper y-wall
        if( y(i) .gt. (boxy/2.d0 - radius(i)) ) then
          imagex = x(i)
          imagey = boxy - y(i)

          !__________________________________
          rx = x(i) - imagex
          ry = y(i) - imagey

          r2 = (rx**2 + ry**2)

          ! assume the particle and its image have the same size
          sigma = radius(i) + radius(i)

          r = dsqrt(r2)

          ff = (1.d0 - r/sigma)

          wallfx(i) = wallfx(i) + (rx/r)*ff/sigma
          wallfy(i) = wallfy(i) + (ry/r)*ff/sigma

          wallV = wallV + 0.5d0*ff**2
          !__________________________________
        endif

      enddo !do i=1,npart

      wallV = wallV/dble(npart)  ! wallV is wall P.E. per particle

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine inclusion_force(x,y,i_x,i_y,fx,fy,i_fx,i_fy,npart,n_i,
     &                           boxx,boxy,radius,i_radius,
     &                           intake_range,radius1)

      implicit none
      integer i,j,k
      integer npart,maxnpart
      parameter(maxnpart=4096)
      double precision x(maxnpart),y(maxnpart),fx(maxnpart),fy(maxnpart)
      double precision boxx,boxy
      double precision radius(maxnpart)
      double precision rx,ry,r,r2,sigma,sigma2,ff
      integer n_i,maxn_i ! n_i = num of inclusion
      parameter(maxn_i=10)
      double precision i_x(maxn_i),i_y(maxn_i)
      double precision i_fx(maxn_i),i_fy(maxn_i)
      double precision i_radius(maxn_i)
      double precision intake_range,radius1,acting_i_radius
      integer dir_switch !0: right; 1: up; 2: down
      common /motion_mode/ dir_switch
      ! L: left; R: right; U: up; D: down walls
      double precision inlet_wallL,inlet_wallR,inlet_wallU,inlet_wallD

      do i=1,n_i
        i_fx(i) = 0.d0
        i_fy(i) = 0.d0
      enddo

      do i=1,npart
        fx(i) = 0.d0
        fy(i) = 0.d0
      enddo

      do i=1,n_i

        if(dir_switch .eq. 0) then     ! roving right
          inlet_wallL = i_x(1) + i_radius(1)-intake_range*radius1
          inlet_wallR = i_x(1) + i_radius(1)
          inlet_wallD = i_y(1) - intake_range/2.d0*radius1
          inlet_wallU = i_y(1) + intake_range/2.d0*radius1
        elseif(dir_switch .eq. 1) then ! roving up
          inlet_wallL = i_x(1) - intake_range/2.d0*radius1
          inlet_wallR = i_x(1) + intake_range/2.d0*radius1
          inlet_wallD = i_y(1) + i_radius(1)-intake_range*radius1
          inlet_wallU = i_y(1) + i_radius(1)
        elseif(dir_switch .eq. 2) then ! roving down
          inlet_wallL = i_x(1) - intake_range/2.d0*radius1
          inlet_wallR = i_x(1) + intake_range/2.d0*radius1
          inlet_wallD = i_y(1) - i_radius(1)
          inlet_wallU = i_y(1) - i_radius(1)+intake_range*radius1
        endif

        do j=1,npart

          if(dir_switch .eq. 0) then     ! roving right
            if(     x(j) .ge. inlet_wallL
     &        .and. y(j) .ge. inlet_wallD
     &        .and. y(j) .le. inlet_wallU) then
              acting_i_radius = i_radius(1) - intake_range*radius1
            else
              acting_i_radius = i_radius(1)
            endif
          elseif(dir_switch .eq. 1) then ! roving up
            if(     y(j) .ge. inlet_wallD
     &        .and. x(j) .ge. inlet_wallL
     &        .and. x(j) .le. inlet_wallR) then
              acting_i_radius = i_radius(1) - intake_range*radius1
            else
              acting_i_radius = i_radius(1)
            endif
          elseif(dir_switch .eq. 2) then ! roving down
            if(     y(j) .le. inlet_wallU
     &        .and. x(j) .ge. inlet_wallL
     &        .and. x(j) .le. inlet_wallR) then
              acting_i_radius = i_radius(1) - intake_range*radius1
            else
              acting_i_radius = i_radius(1)
            endif
          endif

          rx = i_x(i) - x(j)
          ry = i_y(i) - y(j)
          r2 = (rx*rx + ry*ry)
          
          sigma = acting_i_radius + radius(j)
          sigma2 = sigma * sigma

          if(r2 .lt. sigma2) then

            r = dsqrt(r2)
           
            ff = (1.d0 - r/sigma)

            i_fx(i) = i_fx(i) + (rx/r)*ff/sigma
            i_fy(i) = i_fy(i) + (ry/r)*ff/sigma
            fx(j) = fx(j) - (rx/r)*ff/sigma
            fy(j) = fy(j) - (ry/r)*ff/sigma

          endif

        enddo !do j=1,npart
      enddo !do i=1,n_i

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine inlet_wallforce(x,y,wallfx,wallfy,wallV,
     &                           i_x,i_y,i_wallfx,i_wallfy,
     &                           radius,npart,i_radius,n_i,
     &                           intake_range,radius1)

      implicit none
      integer i,j,k
      integer npart,maxnpart
      parameter(maxnpart=4096)
      double precision x(maxnpart),y(maxnpart)
      double precision wallfx(maxnpart),wallfy(maxnpart)
      double precision imagex,imagey
      double precision wallV
      double precision radius(maxnpart)
      double precision rx,ry,r,r2,sigma,sigma2,ff
      integer n_i,maxn_i ! n_i = num of inclusion
      parameter(maxn_i=10)
      double precision i_x(maxn_i),i_y(maxn_i)
      double precision i_wallfx(maxn_i),i_wallfy(maxn_i)
      double precision i_radius(maxn_i)
      double precision intake_range,radius1
      integer dir_switch !0: right; 1: up; 2: down
      common /motion_mode/ dir_switch
      ! L: left; R: right; U: up; D: down walls
      double precision inlet_wallL,inlet_wallR,inlet_wallU,inlet_wallD

      wallV = 0.d0

      do i=1,n_i
        i_wallfx(i) = 0.d0
        i_wallfy(i) = 0.d0
      enddo

      do i=1,npart
        wallfx(i) = 0.d0
        wallfy(i) = 0.d0
      enddo

      do i=1,n_i

        if(dir_switch .eq. 0) then     ! roving right
          inlet_wallL = i_x(1) + i_radius(1)-intake_range*radius1
          inlet_wallR = i_x(1) + i_radius(1)
          inlet_wallD = i_y(1) - intake_range/2.d0*radius1
          inlet_wallU = i_y(1) + intake_range/2.d0*radius1
        elseif(dir_switch .eq. 1) then ! roving up
          inlet_wallL = i_x(1) - intake_range/2.d0*radius1
          inlet_wallR = i_x(1) + intake_range/2.d0*radius1
          inlet_wallD = i_y(1) + i_radius(1)-intake_range*radius1
          inlet_wallU = i_y(1) + i_radius(1)
        elseif(dir_switch .eq. 2) then ! roving down
          inlet_wallL = i_x(1) - intake_range/2.d0*radius1
          inlet_wallR = i_x(1) + intake_range/2.d0*radius1
          inlet_wallD = i_y(1) - i_radius(1)
          inlet_wallU = i_y(1) - i_radius(1)+intake_range*radius1
        endif

        do j=1,npart

          if(     x(j) .ge. inlet_wallL
     &      .and. x(j) .le. inlet_wallR
     &      .and. y(j) .ge. inlet_wallD
     &      .and. y(j) .le. inlet_wallU) then

            !if ! NOT roving left
              ! hit left (L) inlet-wall
              if( x(j) .lt. (inlet_wallL + radius(j)) ) then
                imagex = inlet_wallL - (x(j) - inlet_wallL)
                imagey = y(j)

                !__________________________________
                rx = x(j) - imagex
                ry = y(j) - imagey

                r2 = (rx**2 + ry**2)

                ! assume the particle and its image have the same size
                sigma = radius(j) + radius(j)

                r = dsqrt(r2)

                ff = (1.d0 - r/sigma)

                wallfx(j) = wallfx(j) + (rx/r)*ff/sigma
                wallfy(j) = wallfy(j) + (ry/r)*ff/sigma
                i_wallfx(i) = i_wallfx(i) - (rx/r)*ff/sigma
                i_wallfy(i) = i_wallfy(i) - (ry/r)*ff/sigma

                wallV = wallV + 0.5d0*ff**2
                !__________________________________
              endif
            !endif

            if(dir_switch .ne. 0) then ! NOT roving right
              ! hit right (R) inlet-wall
              if( x(j) .gt. (inlet_wallR - radius(j)) ) then
                imagex = inlet_wallR + (inlet_wallR - x(j))
                imagey = y(j)

                !__________________________________
                rx = x(j) - imagex
                ry = y(j) - imagey

                r2 = (rx**2 + ry**2)

                ! assume the particle and its image have the same size
                sigma = radius(j) + radius(j)

                r = dsqrt(r2)

                ff = (1.d0 - r/sigma)

                wallfx(j) = wallfx(j) + (rx/r)*ff/sigma
                wallfy(j) = wallfy(j) + (ry/r)*ff/sigma
                i_wallfx(i) = i_wallfx(i) - (rx/r)*ff/sigma
                i_wallfy(i) = i_wallfy(i) - (ry/r)*ff/sigma

                wallV = wallV + 0.5d0*ff**2
                !__________________________________
              endif
            endif

            if(dir_switch .ne. 2) then ! NOT roving down
              ! hit lower/down (D) inlet-wall
              if( y(j) .lt. (inlet_wallD + radius(j)) ) then
                imagex = x(j)
                imagey = inlet_wallD - (y(j) - inlet_wallD)

                !__________________________________
                rx = x(j) - imagex
                ry = y(j) - imagey

                r2 = (rx**2 + ry**2)

                ! assume the particle and its image have the same size
                sigma = radius(j) + radius(j)

                r = dsqrt(r2)

                ff = (1.d0 - r/sigma)

                wallfx(j) = wallfx(j) + (rx/r)*ff/sigma
                wallfy(j) = wallfy(j) + (ry/r)*ff/sigma
                i_wallfx(i) = i_wallfx(i) - (rx/r)*ff/sigma
                i_wallfy(i) = i_wallfy(i) - (ry/r)*ff/sigma

                wallV = wallV + 0.5d0*ff**2
                !__________________________________
              endif
            endif

            if(dir_switch .ne. 1) then ! NOT roving up
              ! hit upper/up (U) inlet-wall
              if( y(j) .gt. (inlet_wallU - radius(j)) ) then
                imagex = x(j)
                imagey = inlet_wallU + (inlet_wallU - y(j))

                !__________________________________
                rx = x(j) - imagex
                ry = y(j) - imagey

                r2 = (rx**2 + ry**2)

                ! assume the particle and its image have the same size
                sigma = radius(j) + radius(j)

                r = dsqrt(r2)

                ff = (1.d0 - r/sigma)

                wallfx(j) = wallfx(j) + (rx/r)*ff/sigma
                wallfy(j) = wallfy(j) + (ry/r)*ff/sigma
                i_wallfx(i) = i_wallfx(i) - (rx/r)*ff/sigma
                i_wallfy(i) = i_wallfy(i) - (ry/r)*ff/sigma

                wallV = wallV + 0.5d0*ff**2
                !__________________________________
              endif
            endif

          endif

        enddo !do j=1,npart
      enddo !do i=1,n_i

      wallV = wallV/dble(npart)  ! wallV is wall P.E. per particle

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine self_propel_force(vx,vy,self_propel_fx,self_propel_fy,
     &                             npart)

      implicit none
      integer i,j,k
      integer npart,maxnpart
      parameter(maxnpart=4096)
      double precision vx(maxnpart),vy(maxnpart)
      double precision v(maxnpart),target_v0
      double precision mu
      common /self_propel/ mu
      double precision self_propel_fx(maxnpart),self_propel_fy(maxnpart)

      target_v0 = 0.d0

      do i=1,npart
        v(i) = dsqrt(vx(i)*vx(i) + vy(i)*vy(i))

        if(v(i) .gt. 0.d0) then
          self_propel_fx(i) = mu * (target_v0 - v(i)) * vx(i) / v(i)
          self_propel_fy(i) = mu * (target_v0 - v(i)) * vy(i) / v(i)
        else
          self_propel_fx(i) = 0.d0
          self_propel_fy(i) = 0.d0
        endif
      enddo
    
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine inclusion_self_propel_force(i_vx,i_vy,
     &                                       i_self_propel_fx,
     &                                       i_self_propel_fy,n_i)

      implicit none
      integer i,j,k
      integer n_i,maxn_i ! n_i = num of inclusion
      parameter(maxn_i=10)
      double precision i_vx(maxn_i),i_vy(maxn_i)
      double precision i_v(maxn_i),target_v0
      double precision mu
      common /self_propel/ mu
      double precision i_self_propel_fx(maxn_i),i_self_propel_fy(maxn_i)

      target_v0 = 0.d0

      do i=1,n_i
        i_v(i) = dsqrt(i_vx(i)*i_vx(i) + i_vy(i)*i_vy(i))

        if(i_v(i) .gt. 0.d0) then
          i_self_propel_fx(i) = mu*(target_v0 - i_v(i))* i_vx(i)/ i_v(i)
          i_self_propel_fy(i) = mu*(target_v0 - i_v(i))* i_vy(i)/ i_v(i)
        else
          i_self_propel_fx(i) = 0.d0
          i_self_propel_fy(i) = 0.d0
        endif
      enddo
    
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine KE_cal(vx,vy,mass,i_vx,i_vy,i_mass,tot_ke,npart,n_i)

      implicit none
      integer i,j,k
      integer npart,maxnpart
      parameter(maxnpart=4096)
      double precision vx(maxnpart),vy(maxnpart),mass(maxnpart)
      double precision ke
      integer n_i,maxn_i ! n_i = num of inclusion
      parameter(maxn_i=10)
      double precision i_vx(maxn_i),i_vy(maxn_i),i_mass(maxn_i)
      double precision i_ke
      double precision tot_ke

      ! calculate K.E. per particle
      ke = 0.d0

      do i=1,npart
        ke = ke + 0.5d0 * mass(i) * (vx(i)*vx(i) + vy(i)*vy(i))
      enddo

      ke = ke / dble(npart)  ! ke is K.E. per particle

      ! calculate K.E. per inclusion
      i_ke = 0.d0

      do i=1,n_i
        i_ke = i_ke
     &       + 0.5d0 * i_mass(i) * (i_vx(i)*i_vx(i) + i_vy(i)*i_vy(i))
      enddo

      i_ke = i_ke / dble(n_i)  ! i_ke is K.E. per inclusion

      ! calculate K.E. per (particle + inclusion)
      tot_ke = ke + i_ke

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      double precision ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION gasdev(idum)
      INTEGER idum
      double precision gasdev
CU    USES ran2
      INTEGER iset
      double precision fac,gset,rsq,v1,v2,ran2
      SAVE iset,gset
      DATA iset/0/
      if (idum.lt.0) iset=0
      if (iset.eq.0) then
1       v1=2.*ran2(idum)-1.
        v2=2.*ran2(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
