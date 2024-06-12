import numpy as np

def vis_mesh(sim,nx,nz,dn,lovera):
        
    x = np.linspace(0,2*np.pi/dn,nx)    # r
    z = np.linspace(0,2*np.pi,nz)-np.pi # theta

    xp = np.zeros([nx,nz]) # R
    yp = np.zeros([nx,nz]) # Z 
    zp = np.zeros([nx,nz]) # phi=0 plane

    # minor radius
    r = sim.rmin+(dn*x/(2*np.pi)-0.5)*lovera

    # MXH angle
    a = z + sim.shape_cos[0]
    + sim.shape_cos[1]*np.cos(z) 
    + sim.shape_cos[2]*np.cos(2*z) 
    + sim.shape_cos[3]*np.cos(3*z)
    + sim.shape_cos[4]*np.cos(4*z)
    + sim.shape_cos[5]*np.cos(5*z)
    + sim.shape_cos[6]*np.cos(6*z)
    + np.arcsin(sim.delta)*np.sin(z) 
    - sim.zeta*np.sin(2*z) 
    + sim.shape_sin[3]*np.sin(3*z)
    + sim.shape_sin[4]*np.sin(4*z)
    + sim.shape_sin[5]*np.sin(5*z)
    + sim.shape_sin[6]*np.sin(6*z)
    
    # MESH
    xp[:,:] = sim.rmaj+r[:,None]*np.cos(a[None,:])
    yp[:,:] = sim.zmag+sim.kappa*r[:,None]*np.sin(z[None,:])

    return x,z,xp,yp,zp

def vis_geo(sim,geo):
    
    geo.signb_in=1 # fix
    geo.geo_rmin_in=sim.rmin
    geo.geo_rmaj_in=sim.rmaj
    geo.geo_drmaj_in=sim.shift
    geo.geo_zmag_in=sim.zmag
    geo.geo_dzmag_in=sim.dzmag
    geo.geo_q_in=sim.q
    geo.geo_s_in=sim.shear
    geo.geo_kappa_in=sim.kappa
    geo.geo_s_kappa_in=sim.s_kappa
    geo.geo_beta_star_in=sim.beta_star
    # Antisymmetric
    geo.geo_shape_cos0_in = sim.shape_cos[0]
    geo.geo_shape_cos1_in = sim.shape_cos[1]
    geo.geo_shape_cos2_in = sim.shape_cos[2]
    geo.geo_shape_cos3_in = sim.shape_cos[3]
    geo.geo_shape_cos4_in = sim.shape_cos[4]
    geo.geo_shape_cos5_in = sim.shape_cos[5]
    geo.geo_shape_cos6_in = sim.shape_cos[6]
    # Symmetric
    geo.geo_delta_in=sim.delta
    geo.geo_zeta_in=sim.zeta
    geo.geo_shape_sin3_in = sim.shape_sin[3]
    geo.geo_shape_sin4_in = sim.shape_sin[4]
    geo.geo_shape_sin5_in = sim.shape_sin[5]
    geo.geo_shape_sin6_in = sim.shape_sin[6]
    # Derivative of antisymmetric
    geo.geo_shape_s_cos0_in = sim.shape_s_cos[0]
    geo.geo_shape_s_cos1_in = sim.shape_s_cos[1]
    geo.geo_shape_s_cos2_in = sim.shape_s_cos[2]
    geo.geo_shape_s_cos3_in = sim.shape_s_cos[3]
    geo.geo_shape_s_cos4_in = sim.shape_s_cos[4]
    geo.geo_shape_s_cos5_in = sim.shape_s_cos[5]
    geo.geo_shape_s_cos6_in = sim.shape_s_cos[6]
    # Derivative of symmetric
    geo.geo_s_delta_in=sim.s_delta
    geo.geo_s_zeta_in=sim.s_zeta
    geo.geo_shape_s_sin3_in = sim.shape_s_sin[3]
    geo.geo_shape_s_sin4_in = sim.shape_s_sin[4]
    geo.geo_shape_s_sin5_in = sim.shape_s_sin[5]
    geo.geo_shape_s_sin6_in = sim.shape_s_sin[6]
    
    return geo
