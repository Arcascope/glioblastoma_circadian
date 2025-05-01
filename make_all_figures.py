import pickle

import utils
from fig_compliance import make_compliance_figure
from fig_histograms import make_histograms
from fig_scatter_plot import make_scatter_plot_figure, combine_actogram_and_scatterplot
from fig_timeline import make_timeline_figure
from do_statistics import do_statistics
from make_grids import make_grids

if __name__ == '__main__':
    with open(utils.snapshot_path, 'rb') as input_file:
        people = pickle.load(input_file)
    do_statistics(people)
    make_timeline_figure()
    make_histograms(people)
    make_compliance_figure(people)
    make_scatter_plot_figure(people)
    combine_actogram_and_scatterplot()

    for person in people:
        make_scatter_plot_figure([person], save_name=person.subject_id)
    make_grids()
