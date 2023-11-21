import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import ScalarMappable
import pandas as pd

'''
A script which takes TP output from TPStreamer and plots the found TPs to a pdf file.
Example usage : 
    python ThresholdScan.py -i TP_AbsRS_t25_protoDUNEnoise_n10.txt -e 6 -t 25 30 40 50 -o ThresholdScan_AbsRS_protoDUNEnoise.pdf

Possible to do a threshold scan - i.e. see what the TPs would look like for various thresholds (strictly higher than the ones the TPs were generated with in LArSoft). 
The TPs are plotted in channel-tick parameter space, with the size of the marker corresponding to hit TOT (longer hits are bigger)
and the colour corresponding to the hit SADC (summed ADC). The hits are overlaid on top of tri-colored squares which (hopefully) 
correspond to the MC producer for that hit: 
            white squares == noise hits,
            red squares == signal hits (SN neutrino hits) 
            blue squares  == background hits (hits produced by radiological backgrounds - actual physics signal, albeit undesirable) 

Side note: the goal of the TPGen algorithms is to filter physics signal hits (neutrino + radiological) from the noise hits.
The goal of Trigger algorithms is to further filter actual signal hits (neutrino only) from backgrond hits (radiological). 
'''

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Plot TP data with various thresholds.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input file")
    parser.add_argument("-e", "--event", type=int, required=True, help="Event to plot")
    parser.add_argument("-t", "--thresholds", type=int, nargs='+', required=True, help="Thresholds list (space-separated)")
    parser.add_argument("-o", "--output", type=str, help="Output filename")

    args = parser.parse_args()

    # Load the hit data and put it into a pandas dataframe
    hits_ = pd.read_csv(args.input, delimiter=',', usecols=range(10), names=['event', 'plane', 'startT', 'endT', 'peakT', 'TOT', 'channel', 'SADC', 'peakADC', 'ID'])

    # Colormaps and plot settings
    plt.rcParams.update({'font.size': 7})
    plt.rcParams["axes.facecolor"] = 'gainsboro'

    cmap = plt.cm.colors.ListedColormap(['white', 'red', 'blue'])
    cmap_hits = 'seismic'

    # Define hit tagging marker information
    markers_info = {
        'noise': {'marker': 's', 'color': 'white', 'position': (0.1, 0.95)},
        'signal': {'marker': 's', 'color': 'red', 'position': (0.1, 0.9)},
        'radio. bgd.': {'marker': 's', 'color': 'blue', 'position': (0.1, 0.85)}
    }

    # PDF file for saving the figures
    if args.output:
        pdf_filename = args.output
    else:
        pdf_filename = "ThresholdScan.pdf"

    pdf_pages = PdfPages(pdf_filename)
        
    #Loop over the thresholds and plot the TPs 
    for t in args.thresholds:
        fig, axs = plt.subplots(3, 1, figsize=(6, 8), dpi=200)
        label = f'Threshold = {t} ADC'
        labels = [
            f'Induction plane (U) | {label}',
            f'Induction plane (V)  | {label}',
            f'Collection plane (X)  | {label}'
        ]

        for i in range(3):
            d = hits_[(hits_.event == args.event) & (hits_.plane == i) & (hits_.peakADC > t)]
            # Plot the tagger
            axs[i].scatter(d.peakT, d.channel, s=d.TOT * 1.5,
                           c=d.ID, # Marker color denotes producer ID
                           cmap=cmap, vmin=0, vmax=2,
                           marker='s', alpha=0.5)
            # Plot the hits
            axs[i].scatter(d.peakT, d.channel, s= d.TOT,# Position of marker = [tick, channel], marker size is proportional to TOT
                           c=d.SADC,  # Marker color denotes hit SADC
                           cmap=cmap_hits,
                           marker='o',
                           edgecolor='white',
                           linewidth=0.4)
            axs[i].set_title(labels[i])
            axs[i].set_xlabel("time [tick]")
            axs[i].set_ylabel("channel")
            axs[i].grid(color='white', linestyle='dashed', alpha=0.3)

            # Colorbar
            sm = ScalarMappable(cmap=cmap_hits, norm=plt.Normalize(vmin=d.SADC.min(), vmax=d.SADC.max()))
            sm.set_array([])
            colorbar = plt.colorbar(sm, ax=axs[i])
            colorbar.set_label('Hit SADC')

            # Hit Tagger key
            for label, info in markers_info.items():
                axs[i].scatter([], [], marker=info['marker'], color=info['color'], label=label)
            axs[i].legend(loc='lower center', ncol=3, facecolor='darkgray', labelcolor='k')

            # Center the plot around signal hits
            d = hits_[(hits_.event == args.event) & (hits_.ID == 1)]
            t_offset = 100;  ch_offset = 6000
            axs[i].set_xlim(d.peakT.min() - t_offset, d.peakT.max() + t_offset)
            axs[i].set_ylim(d.channel.min() - ch_offset, d.channel.max() + ch_offset)

        plt.tight_layout()
        pdf_pages.savefig(fig)  

    pdf_pages.close()

    # Show the saved figures (optional)
    # Don't uncomment if plotting for many thresholds as the plots will open in separate windows!
    #plt.show()

if __name__ == "__main__":
    main()
