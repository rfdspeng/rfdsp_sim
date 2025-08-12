import math

def modulation_mapper(bitstream,modorder):
    if modorder == 4:
        tones = 1-2*bitstream[0::2] + 1j*(1-2*bitstream[1::2])
        tones = tones/math.sqrt(2)
    elif modorder == 16:
        tones = (1-2*bitstream[0::4]) * (2-(1-2*bitstream[2::4])) + 1j * (1-2*bitstream[1::4]) * (2-(1-2*bitstream[3::4]))
        tones = tones/math.sqrt(10)
    elif modorder == 64:
        tones = (1-2*bitstream[0::6]) * (4-(1-2*bitstream[2::6])*(2-(1-2*bitstream[4::6]))) + 1j * (1-2*bitstream[1::6]) * (4-(1-2*bitstream[3::6])*(2-(1-2*bitstream[5::6])))
        tones = tones/math.sqrt(42)
    elif modorder == 256:
        tones = (1-2*bitstream[0::8]) * (8-(1-2*bitstream[2::8])*(4-(1-2*bitstream[4::8])*(2-(1-2*bitstream[6::8])))) + 1j * (1-2*bitstream[1::8]) * (8-(1-2*bitstream[3::8])*(4-(1-2*bitstream[5::8])*(2-(1-2*bitstream[7::8]))))
        tones = tones/math.sqrt(170)
    else:
        raise ValueError("modorder must be in [4, 16, 64, 256]")
        
    return tones