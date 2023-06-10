use std::iter::zip;


fn accn(p: &[[f64;3]], gm: &[f64], a: &mut [[f64;3]]) -> () {
    
    for i in a.iter_mut() {
        for j in i.iter_mut() { 
            *j = 0.0;
        }
    }
    //I think the performance detriment is in the inefficient zip/for_eaching in this loop, but
    //need to profile
    for (i,p1) in p.iter().enumerate() {
        for (j,p2) in p[i+1..].iter().enumerate() {
 
            /*let d_ = zip(p1,p2).map( |(x,y)| *x-*y );
            let n: f64 = d_.clone().map( |x| x*x).sum();
            let d_:Vec<_> = d_.map( |x| x/(n*n.sqrt())).collect(); */
            let mut d_ :[f64;3] = [ p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2] ];
            let n = d_[0]*d_[0] + d_[1]*d_[1] + d_[2]*d_[2];
            for dd in d_.iter_mut() {
                *dd = *dd/(n*n.sqrt());
            }
            a[j+i+1].iter_mut().zip(d_.iter()).for_each( |(x,y)| *x += y*gm[i]);
            a[i].iter_mut().zip(d_.iter()).for_each(|(x,y)| *x -= gm[j+i+1]*y);
        }           
    }
}

const w0:f64 = -1.702414383919315268;
const w1:f64 = 1.351207191959657634;
const c:[f64;4] = [w1/2.0, (w0+w1)/2.0, (w0+w1)/2.0, w1/2.0];
const d:[f64;3] = [w1, w0, w1];
const GSCALE:f64 = 2.22972472095e-15;

fn integrate(p: &mut [[f64;3]], v: &mut [[f64;3]], gm: &[f64], a: &mut [[f64;3]], mut accn: impl FnMut(&[[f64;3]],&[f64], &mut [[f64;3]])-> (), dt: f64) -> () {
    for (cc,dd) in zip(c,d)  {
        p.iter_mut().zip(v.iter()).for_each( |(x,y)| zip(x,y).for_each(|(x,y)| *x += *y*cc*dt));
        accn(p, gm, a);
        v.iter_mut().zip(a.iter()).for_each( |(x,y)| zip(x,y).for_each(|(x,y)| *x += *y*dd*dt));
    }
    p.iter_mut().zip(v.iter()).for_each( |(x,y)| zip(x,y).for_each(|(x,y)| *x += *y*c[3]*dt));    
}


fn main() {
    let mut p: [[f64;3]; 11] = [  [2.950832139557744e-3,  -5.425470573959765e-3,  -7.383386124694714e-5],/*Sol*/
        [-3.332497263005759,4.112363636591643,5.891307575840340e-2],/*Jupiter*/
        [-6.535507530500684,6.380917893848125,1.428260570382428e-1],/*Saturn*/
        [2.510983241619308e-1,1.950982735843640e-1,-6.755395306466706e-3],/*Mercury*/
        [-2.620890874732100e-1,-6.808075673769205e-1,6.655856110641488e-3],/*Venus*/
        [9.630524101520794e-1,-3.128557379506953e-1,-2.422339413846039e-4],/*Earth*/
        [9.612598229631927e-1,-3.112560037549189e-1,-2.786213887250675e-4],/*Moon*/
        [-4.415739460116890e-1,1.541057808775840,4.336870492890024e-2],/*Mars*/
        [1.718669476007500e1,9.933497426989296,-1.866007100527533e-1],/*Uranus*/
        [2.615297960363933e1,-1.469969306823841e1,-2.995311280685409e-1],/*Neptune*/
        [6.391931737411634e-1,-1.121186352054070e-1,1.870482953942138e-1]   ];
    let mut v: [[f64;3]; 11] = [   [	6.838345177814781e-6,5.026348301755031e-6,-2.071993311542051e-7],/*Sol*/
        [-5.946882707429968e-3,-4.397581117062436e-3,1.511463721442163e-4],/*Jupiter*/
        [-4.208760152990489e-3,-4.000694588560874e-3,2.380878809848479e-4],/*Saturn*/
        [-2.319387349809163e-2,2.314935923144703e-2,4.029155127823395e-3],/*Mercury*/
        [1.869504930197275e-2,-7.474624825638375e-3,-1.179108325817677e-3],/*Venus*/
        [4.975052051056192e-3,1.633653719659161e-2,1.120657612318169e-5],/*Earth*/
        [4.561144105343879e-3,1.586231559918362e-2,6.830691822388147e-5],/*Moon*/
        [-1.288189187342671e-2,-2.690962310018411e-3,2.691923720938971e-4],/*Mars*/
        [-2.000317996107281e-3,3.225767036414139e-3,3.845761974599044e-5],/*Uranus*/
		[1.517663638954640e-3,2.753873265004750e-3,-9.168067707913533e-5],/*Neptune*/
        [-1.551860810762931e-2,-2.496616444169774e-2,4.202270253608422e-4]  ];
    let mut a : [[f64;3]; 11] = [[0.0;3];11];
    
    let gm = [GSCALE*132712440041.93938,GSCALE*126686531.900,GSCALE*37931206.234,GSCALE*22031.86855, GSCALE*324858.592, GSCALE*398600.435436, GSCALE*4902.800066 /*Moon*/,
				    GSCALE*42828.375214, GSCALE*5793951.256, GSCALE*6835099.97, 0.0];

    let mut t = 0.0;
    let dt = 0.1;
    let max_t = 365.25*200.0;
    let mut dist:f64 = zip(p[10],p[0]).map(|(x,y)| (x-y)*(x-y)).sum();
    let mut olddist = dist;
    let mut oldolddist = dist;


    while t<max_t {
        oldolddist = olddist;
        olddist = dist;

        integrate(&mut p, &mut v, &gm, &mut a, &accn, dt);

        dist = zip(p[10],p[0]).map(|(x,y)| (x-y)*(x-y)).sum();


        if ( oldolddist > olddist) && (dist > olddist) {
            println!("Perihelion at: {t}");
        }

        t+=dt;
    }

}
