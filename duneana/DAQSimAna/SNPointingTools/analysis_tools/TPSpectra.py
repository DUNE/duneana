import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["axes.facecolor"]='gainsboro'
plt.rcParams.update({'font.size': 10})


'''
Simple script which takes TP output from TPStreamer and plots the ID-tagged hit spectra. 

Example usage:
python TPSpectra.py -i TP_standardHF_t20_protoDUNEnoise_n10.txt -e -1 -o TPSpectra_standardHF_protoDUNEnoise.pdf


'''

def main():
    parser = argparse.ArgumentParser(description="TP properties histograms from a tpstream file.")
    parser.add_argument("-i", "--input", required=True, help="Input data file")
    parser.add_argument("-e", "--event", type=str, default="-1", help="Event or range of events to plot (use -1 to plot all events)")
    parser.add_argument("-t", "--threshold", help='Specify whether only hit spectra above a certain ADC threshold should be plotted')
    parser.add_argument("-o", "--output", default="TPSpectra.pdf", help="Output PDF filename (optional)")
    
    args = parser.parse_args()

    # Load data
    hits_ = pd.read_csv(args.input, delimiter=',', usecols=range(10), 
                        names=['event','plane','startT','endT','peakT','TOT','channel', 'SADC','peakADC','ID'])

    # Histogram labels
    hits_['ID'] = hits_['ID'].replace(0, 'noise')
    hits_['ID'] = hits_['ID'].replace(1, 'signal')
    hits_['ID'] = hits_['ID'].replace(2, 'radio. bgd')
    
    custom_palette = {
    'noise': 'white',
    'signal': 'red',
    'radio. bgd': 'blue'
    }


    # Filter events based on passed flags
    if args.event != "-1":
        # Convert event argument to a list of event numbers
        events_to_plot = [int(e) for e in args.event.split(',')]
        hits_ = hits_[hits_['event'].isin(events_to_plot)]

    if args.threshold:
        hits_ = hits_[(hits_.peakADC > int(args.threshold))]

    # Define the columns to be plotted and create subplots
    cols = ['SADC', 'peakADC', 'TOT', 'channel', 'peakT', 'plane']
    n_subplots = len(cols)
    n_rows = (n_subplots + 1) // 2
    fig, axes = plt.subplots(n_rows, 2, figsize=(10, 4.5 * n_rows))  
    axes = axes.flatten()

    # Make sure the smallest histograms are always in front
    id_order = hits_.groupby('ID').size().sort_values().index

    for i, column in enumerate(cols):
        ax = axes[i]
        if column == 'plane':
            bin_edges = [-0.4, 0.4, 0.6, 1.4, 1.6, 2.4]
            sns.histplot(data=hits_, x=column, hue='ID', element='step', hue_order=id_order,
                         ax=ax, bins=bin_edges, common_norm=True, palette=custom_palette)
        else:
            sns.histplot(data=hits_, x=column, hue='ID', element='step', hue_order=id_order,
                         ax=ax, bins=60, common_norm=True, palette=custom_palette)
        ax.set_title(f'{column}')
        ax.set_ylabel('Frequency')
        ax.set_yscale('log')

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, format='pdf')
    plt.show()

if __name__ == "__main__":
    main()
