import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

c = 299792458

class lensClass():
    def __init__(self, focLen, pos=None):
        self.focLen = focLen
        self.pos = pos
        

def main():
    freqs = np.linspace(0.3, 1, 5) * 1e12
    lambs = c / freqs

    w0 = 3e-3 #initial beam waist (m)

    #lens setup (m)
    lenses = []
    lenses.append(lensClass(50e-3, 50e-3))
    lenses.append(lensClass(100e-3))

    #For the 2 lens example:
    lens2OptimalPosition = lenses[0].focLen + lenses[1].focLen + lenses[0].pos
    lenses[1].pos = lens2OptimalPosition

    scanLength = lenses[-1].pos + lenses[-1].focLen * 2

    pos = []
    posWidth = []
    for freqInd, freq in enumerate(freqs):
        posT, posWidthT = gaussian_propagation(freq, w0, 0, lenses, scanLength)
        pos.append(np.array(posT))
        posWidth.append(np.array(posWidthT))
        

    fig, ax1 = plt.subplots(1)
    lineLeg = []
    for freqInd, freq in enumerate(freqs):
        lineLeg.append(ax1.plot(pos[freqInd] * 1e3, posWidth[freqInd] * 1e3, linewidth=0.85, label=str(freq*1e-12)+' THz')[0])
        ax1.plot(pos[freqInd] * 1e3, posWidth[freqInd] * -1e3, linewidth=0.85, c = lineLeg[-1].get_c())

    for lens in lenses:
        ax1.plot(np.array([lens.pos, lens.pos]) * 1e3, np.array([-15, 15]) * 1e3, 'k--', linewidth=0.65)
        ax1.arrow(lens.pos * 1e3, 0, lens.focLen * 1e3, 0, color='k', linewidth=0.65, head_width=0.5, head_length=10, length_includes_head=True, alpha=0.2)
        
    ax1.axis([0, scanLength * 1e3, -15, 15])
    ax1.set_xlabel('Scan Dist (mm)')
    ax1.set_ylabel('Beam Radius (mm)')

    plt.legend(handles = lineLeg, loc = 2, fontsize = 8)
    plt.show()


def gaussian_focusing(freq, waistPre, waistPosPre, lensFoc):
# Gaussian beam focusing by a thin lens (see [1])
# Inputs: frequency 
#         waistPre = waist before the lens;
#         waisPosPre = distance of waist before lens (>0: before lens)
#         lensFoc = focal length of lens (>0: converging lens)
# Output: waistPost = new waist;
#         waistPosPost = position of waist from lens
# [1] http://www.mellesgriot.com/products/optics/gb_2_3.htm
    zR = np.pi * waistPre**2 * freq / c   # original rayleigh range
    waistPosPost = lensFoc * (1 + (waistPosPre / lensFoc - 1) / ((waistPosPre / lensFoc - 1)**2 + (zR / lensFoc)**2));  # new waist location (Eq. (9b))
    waistPost = waistPre / np.sqrt((1 - waistPosPre / lensFoc)**2 + (zR / lensFoc)**2);
    return waistPost, waistPosPost


def gaussian_propagation(freq, beamWaist, waistPos, lensList, totDist):
# propagates a Gaussian beam through a series of lenses.
# Input:
#   freq           input frequency
#   beamWaist      waist of the beam at the waist positions (wp)
#   waistPos       waist position
#   lensList       list of lenses
#   totDist        total distance to calculate

    waist = [beamWaist]
    waistP = [waistPos]

    for lensInd, lens in enumerate(lensList):

        lensList[lensInd].wPre = waist[lensInd]
        lensList[lensInd].wPrePos = lens.pos - waistP[lensInd]

        waistPost, waistPosPost = gaussian_focusing(freq, lens.wPre, lens.wPrePos, lens.focLen)
        
        lensList[lensInd].wPost = waistPost
        lensList[lensInd].wPostPos = lens.pos + waistPosPost
        
        waist.append(lens.wPost)
        waistP.append(lens.wPostPos)
        
        
    propPos = []
    propWidth = []
    
    #Pre lens scan:
    posPoints = np.linspace(waistPos, lensList[0].pos, 100)
    posWidth = beamWaist * np.sqrt(1 + ((posPoints - waistPos) / (np.pi * beamWaist**2 * freq / c))**2)
    propPos.extend(posPoints)
    propWidth.extend(posWidth)
    
    #Scan width between each lens:
    for lensIndPlus1, lens in enumerate(lensList[1:]):
        posPoints = np.linspace(lensList[lensIndPlus1].pos, lens.pos, 100)
        posWidth = lens.wPre * np.sqrt(1 + ((posPoints - (lens.pos - lens.wPrePos)) / (np.pi * lens.wPre**2 * freq / c))**2)
        propPos.extend(posPoints)
        propWidth.extend(posWidth)
        
    #Post lens scan:
    posPoints = np.linspace(lensList[-1].pos, totDist, 100)
    posWidth = lensList[-1].wPost * np.sqrt(1 + ((posPoints - lensList[-1].wPostPos) / (np.pi * lensList[-1].wPost**2 * freq / c))**2)
    propPos.extend(posPoints)
    propWidth.extend(posWidth)

    return propPos, propWidth

main()