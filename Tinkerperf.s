!-------------------------------------------------------------------
 vpcompressd %zmm3, vec_mp_kglobvec1_(,%rcx,4){%k1}
 vpcompressd %zmm4, vec_mp_kbisvec1_(,%rcx,4){%k1}
 vpcompressd %zmm5, vec_vdw_mp_kvvec1_(,%rcx,4){%k1}
 addl        %edx, %r14d
 incl        %edx
 movslq      %edx, %rdx
 vpcompressd %zmm6, -4+vec_vdw_mp_kvlocvec1_(,%rdx,4){%k1}
!-------------------------------------------------------------------
 jle       ..B2.324
 cmpl      $8, %r13d
 jl        ..B2.460
 movl      %r13d, %edx
 xorl      %eax, %eax
 andl      $-8, %edx
 movslq    %edx, %rdx
 vbroadcastsd vdwpot_mp_ghal_(%rip), %zmm17
 vbroadcastsd vdwpot_mp_dhal_(%rip), %zmm16
 vmovups   -816(%rbp), %zmm20
 movl      %r14d, -456(%rbp)
 movq      %rdx, %r14
 movq      %r12, -152(%rbp)
 movq      %rax, %r12
 # MAIN VECTOR TYPE: 64-bits floating point
 vmovups   vec_vdw_mp_rikvec_(,%r13,8), %zmm19           
 vmovups   vec_vdw_mp_rvvec2_(,%r13,8), %zmm18           
 vmulpd    vec_vdw_mp_rik2vec_(,%r13,8), %zmm19, %zmm5   
 vmulpd    %zmm18, %zmm18, %zmm2                         
 vmulpd    %zmm19, %zmm5, %zmm6                          
 vmulpd    %zmm2, %zmm2, %zmm3                           
 vmulpd    %zmm18, %zmm2, %zmm4                          
 vmovupd   %zmm5, vec_vdw_mp_rik3vec_(,%r13,8)           
 vmulpd    %zmm19, %zmm6, %zmm7                          
 vmulpd    %zmm4, %zmm3, %zmm0                           
 vmovupd   %zmm6, vec_vdw_mp_rik4vec_(,%r13,8)           
 vmulpd    %zmm19, %zmm7, %zmm8                          
 vmovupd   %zmm0, vec_vdw_mp_rv7vec_(,%r13,8)            
 vmovupd   %zmm7, vec_vdw_mp_rik5vec_(,%r13,8)           
 vmulpd    %zmm19, %zmm8, %zmm9                          
 vmovupd   %zmm8, vec_vdw_mp_rik6vec_(,%r13,8)           
 vfmadd213pd %zmm9, %zmm17, %zmm0                        
 vmovupd   %zmm9, vec_vdw_mp_rik7vec_(,%r13,8)           
 vmovaps   %zmm20, %zmm1
 call      *__svml_pown8_z0@GOTPCREL(%rip)
 vfmadd231pd %zmm16, %zmm18, %zmm19
 vmovaps   %zmm20, %zmm1
 vmovupd   %zmm0, vec_vdw_mp_invrhovec_(,%r12,8)
 vmovaps   %zmm19, %zmm0
 call      *__svml_pown8_z0@GOTPCREL(%rip)
 vmovupd   %zmm0, vec_vdw_mp_invtmpvec_(,%r12,8)
 addq      $8, %r12
 cmpq      %r14, %r12
 jb        ..B2.322
!-------------------------------------------------------------------
