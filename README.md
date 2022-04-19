# DualThreshold_BeamTraining
This is fast beam training algorithm based on E. Eltayeb's opportunistic beam training. In this work, two thresholds algorithm is proposed. If beam SNR exceeds the first threshold, we use that beam for communication without further searching. Else if SNR exceeds the second threshold(=0.65 * the first threshold), we search narrower range centered on that beam(range *= 1/4 or 1/8) so that we can find the desired beam.

Many of annotation is written in Korean. If you need English or another language version, please contact me: initial@snu.ac.kr
