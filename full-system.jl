#Note, we *need* StaticArrays here to make this approach C performance - the built in dynamic arrays incur a 10x performance penalty even with preallocation etc

using StaticArrays
using LinearAlgebra

p = @MMatrix [  2.950832139557744e-3  -5.425470573959765e-3  -7.383386124694714e-5 ;
                      -3.332497263005759 4.112363636591643 5.891307575840340e-2 ;
                      -6.535507530500684 6.380917893848125 1.428260570382428e-1 ;
                      2.510983241619308e-1 1.950982735843640e-1 -6.755395306466706e-3 ;
                      -2.620890874732100e-1 -6.808075673769205e-1 6.655856110641488e-3 ;
                      9.630524101520794e-1 -3.128557379506953e-1 -2.422339413846039e-4 ;
                      9.612598229631927e-1 -3.112560037549189e-1 -2.786213887250675e-4 ;
                      -4.415739460116890e-1 1.541057808775840 4.336870492890024e-2 ;
                      1.718669476007500e1 9.933497426989296 -1.866007100527533e-1 ;
                      2.615297960363933e1 -1.469969306823841e1 -2.995311280685409e-1 ;
                      6.391931737411634e-1 -1.121186352054070e-1 1.870482953942138e-1 ; ]

v = @MMatrix [ 6.838345177814781e-6 5.026348301755031e-6 -2.071993311542051e-7 ;
                      -5.946882707429968e-3 -4.397581117062436e-3 1.511463721442163e-4 ; 
                      -4.208760152990489e-3 -4.000694588560874e-3 2.380878809848479e-4 ;
                      -2.319387349809163e-2 2.314935923144703e-2 4.029155127823395e-3 ;
                      1.869504930197275e-2 -7.474624825638375e-3 -1.179108325817677e-3 ;
                      4.975052051056192e-3 1.633653719659161e-2 1.120657612318169e-5 ;
                      4.561144105343879e-3 1.586231559918362e-2 6.830691822388147e-5 ;
                      -1.288189187342671e-2 -2.690962310018411e-3 2.691923720938971e-4 ; 
                      -2.000317996107281e-3 3.225767036414139e-3 3.845761974599044e-5 ;
                      1.517663638954640e-3 2.753873265004750e-3 -9.168067707913533e-5 ;
                      -1.551860810762931e-2 -2.496616444169774e-2 4.202270253608422e-4 ; ]

const Gm = @SVector [ i*2.22972472095e-15 for i in [ 132712440041.93938,126686531.900,37931206.234,22031.86855, 324858.592, 398600.435436, 4902.800066, 42828.375214, 5793951.256, 6835099.97,0.0]]

a = MMatrix{11,3}(zeros(Float64, (11,3)))



const w0::Float64 = -1.702414383919315268


const w1::Float64 = 1.351207191959657634


const c = (w1/2, (w0+w1)/2, (w0+w1)/2, w1/2)

const d = (w1, w0, w1)


function accns!(psn,Gm,acn) 
           a .= 0.0
           d = @MVector [0.0,0.0,0.0]
           lst = size(psn)[1]
           for i in 1:lst
               for j in i+1:lst
                   d .= psn[i,:] .- psn[j,:]
                   n::Float64 = d ⋅ d #sum(d .^ 2) better to make this explicit
                   d ./= (n*sqrt(n))
                   acn[j,:] .+= Gm[i]*d
                   acn[i,:] .-= Gm[j]*d
               end
           end
       end

function integrate!(posns, velocities, Gmasses, acns, accn::Function, dt)
                  if !(size(posns) == size(velocities) )
                      throw(ArgumentError("All vectors passed to integrate!() must have equal length"))
                  end
                  for i in eachindex(d)
                      posns .+= velocities .* (c[i] * dt)
                      accn(posns,Gmasses,acns)
                      velocities .+= acns .* (d[i] * dt)
                  end
                  posns .+= velocities .* (c[4] * dt)
              end


function sim(p,v,a,dt,max_t)
           t = 0
           lst = size(p)[1]
           dist = sum( (p[lst,:] .- p[1,:]) .^ 2 )
           olddist = dist
           while t < max_t
               oldolddist = olddist
               olddist = dist
               
               integrate!(p,v,Gm,a,accns!, dt)
               dist = sum( (p[lst,:] .- p[1,:]) .^ 2 )
               
               if (oldolddist > olddist) && (dist > olddist)
                   println("Perihelion at $t")
               end
               
               t += dt
           end
       end


