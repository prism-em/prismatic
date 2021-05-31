import h5py
import numpy as np
from scipy import signal

###########################
######### Classes #########
###########################

# Class for implementing chromatic aberration
class Cc:
    def __init__(self, images,defocus):
        

        # Inputs are an array of defocus and images 
        self.defocus = np.array(defocus)

        #print(self.defocus)
        self.img = np.array(images)
        
        # Number of images
        self.n = len(self.defocus)

        self.std_dev = []
        for i in range(int((self.n)/2)):
            self.std_dev.append(self.defocus[int(self.n/2)+i+1] - self.defocus[int(self.n/2)-i - 1])
            
        self.std_dev = np.array(self.std_dev)/2

        self.initial_defocus = self.defocus[int(self.n/2)]
        self.sigma = self.std_dev[0]
        self.dz = self.sigma*2.355

        self.i1 = 0
        self.i2 = 0
        self.i3 = 0



    # Gaussian function
    def Gaussian(self, x, sigma, mu):
        return (1/np.sqrt(2*np.pi*sigma**2))* np.exp(-(x - mu)**2/(2*sigma**2))
        
    # Calulate the weights from Gaussian
    def calculateWeights(self, defocus, sigma):

        self.pos = np.where(self.defocus == defocus)[0][0]
        h = self.pos - int(self.n/2)
        weights = np.zeros(self.n)

        for i in range(self.n):
            weights[i] = self.Gaussian(self.defocus[i], sigma, defocus)


        for i in range(2*abs(h)):
            if h > 0:
                weights[i] = 0
            else:
                weights[-(i+1)] = 0


        self.weights = weights
        self.norm = np.sum(self.weights)

        if self.norm ==0:
            self.weights = np.zeros(self.n)

        else:
            self.weights = weights/self.norm
            
        plotdef = []
        plotweights= []
        for i in range(self.n):
            if self.weights[i] != 0:
                plotdef.append(self.defocus[i])
                plotweights.append(self.weights[i])

        self.plotdefocus = np.array(plotdef)
        self.plotweights = np.array(plotweights)
        return self.weights


    # Calculates the CC image by multiplying the weights by the images and summing them
    def calculateImage(self, weights):
        img = np.einsum("ijk, i -> jk", self.img, weights)
        img = np.array(img)
        self.CC_img = img
        return img


    # Combines above functions to output CC image
    def computeImage(self, defocus, sigma):
        img = self.calculateImage( self.calculateWeights(defocus, sigma) )
        return img


###########################
######## Utilities ########
###########################
def rad_gaussian(r, sigma):
    return (0.5*np.pi*sigma)*np.exp(-0.5*(r**2.0)/(sigma**2.0))

def rad_cauchy(r, gamma):
		return (gamma)/((x**2) + gamma**2)

def rad_weights(r, p,  mask, fun=rad_gaussian):
    weights = fun(r, p)
    weights *= mask
    weights /= np.sum(weights) #normalize to sum to 1
    return weights

def calc_aberration(coef, q, qtheta, lmb):
    """
    using prismatic aberration convention

    coef is N x 4 np array of 
    ab_0 = M N mag angle
    angle in degrees
    lmb is electron lambda
    """
    assert(len(coef.shape)==2), "coef must be N x 4 np array"
    assert(coef.shape[1]==4), "need 4 coefficients for each variable"
    chi = np.zeros(q.shape, dtype=complex)
    for ab in coef:
        if (ab[0] >= ab[1]) and ( not ( ab[0] + ab[1]) % 2):
            rad = ab[3] * np.pi / 180.0
            cx = ab[2]*np.cos(ab[1]*rad)
            cy = ab[2]*np.sin(ab[1]*rad)

            chi = chi \
                  + cx*((lmb*q)**ab[0])*np.cos(ab[1]*qtheta) \
                  + cy*((lmb*q)**ab[0])*np.sin(ab[1]*qtheta)

    return chi

def calc_defocus_aberration(q, qTheta, lmb, C1):
    """
    wrapper function for general convention above
    uses defocus in C1 convention
    prepares qTheta, q, and dimensionless coef and returns chi
    """
    coef = np.zeros((1,4))
    coef[0,0] = 2
    coef[0,1] = 0
    coef[0,2] = C1*np.pi/lmb

    return calc_aberration(coef, q, qTheta, lmb)

def apply_aberration(psi, chi):
    """
    input psi in realspace
    chi in fourier space
    """
    kpsi = np.fft.fft2(psi) #bring wave fucntion to fourier space
    a_kpsi = kpsi*np.exp(-1j*chi) #apply aberration in fourier space
    a_psi = np.fft.ifft2(a_kpsi) #transorm aberrated psi back to real space
    return a_psi

def get_lambda(fp):
    """
    use h5py to get attribute since metadata a bit broken with py4DSTEM
    """
    tmp_file = h5py.File(fp, "r")
    w_lambda = tmp_file['4DSTEM_simulation']['metadata']['metadata_0']['original']['simulation_parameters'].attrs['lambda']
    tmp_file.close()
    return w_lambda

def get_pixel_size(fp):
    """
    use h5py to get attribute since metadata a bit broken with py4DSTEM
    """
    tmp_file = h5py.File(fp, "r")
    px_x = tmp_file['4DSTEM_simulation']['metadata']['metadata_0']['original']['simulation_parameters'].attrs['eff_pixel_size_x']
    px_y = tmp_file['4DSTEM_simulation']['metadata']['metadata_0']['original']['simulation_parameters'].attrs['eff_pixel_size_x']
    tmp_file.close()
    return px_x, px_y

def get_tilts(fp):
    F_hrtem_tilts = h5py.File(fp, mode='r')
    tilt_x = np.copy(F_hrtem_tilts['4DSTEM_simulation']['data']['datacubes']['HRTEM_virtual']['dim3'])
    tilt_y = np.copy(F_hrtem_tilts['4DSTEM_simulation']['data']['datacubes']['HRTEM_virtual']['dim4'])
    F_hrtem_tilts.close()
    return tilt_x, tilt_y

def get_rdims(fp, datagroup, dataset):
    f = h5py.File(fp, mode='r')
    rx = np.copy(f['4DSTEM_simulation']['data'][datagroup][dataset]['dim1'])
    ry = np.copy(f['4DSTEM_simulation']['data'][datagroup][dataset]['dim2'])
    f.close()
    return rx, ry

def apply_poisson_noise_STEM(data, dose=10, area=False, probe_step=(1.0,1.0), normalize=True):
    """
    if area:
        dose in e/ang^2
    else:
        dose in e/probe position
        
    all STEM outputs have fractional beam intensities
    probe_step should be tuple of probe steps in units of angstroms
    """
    
    if area:
        probe_area = np.prod(probe_step)
        N = dose*probe_area
    else:
        N = dose
        
    if normalize:
        return np.random.poisson(data*N)/N
    else:
        return np.random.poisson(data*N)

def apply_poisson_noise_HRTEM(data, dims, dose=10, normalize=True):
    """
    dose in e/ang^2
    data should be some array of intensities
    HRTEM output has mean intensity of 1
    tilts are only in last dimension
    """
    px_area = np.prod(dims)/np.prod(data.shape[0:2])
    N = dose*px_area # N events per pixel
    
    if normalize:
        return np.random.poisson(data*N)/N
    else:
        return np.random.poisson(data*N)

def source_size(data, rx, ry, sigma):
    """
    source-size effects for STEM simulations
    """
    ryy, rxx = np.meshgrid(ry,rx)
    rxx = rxx-rx[len(rx)//2-1]
    ryy = ryy-ry[len(ry)//2-1]
    rr = np.sqrt(rxx**2.0 + ryy**2.0)
    source_size_kernel = rad_gaussian(rr, kernel_sigma)
    
    return signal.convolve2d(data, source_size_kernel, 'same')/signal.convolve2d(np.ones(np.shape(data)), source_size_kernel, 'same')