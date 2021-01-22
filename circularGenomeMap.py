from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

#set the arrow sigil properties
def _draw_sigil_arrow(
        self, bottom, center, top, startangle, endangle, strand, **kwargs
    ):
        """Draw ARROW sigil (PRIVATE)."""
        if strand == 1:
            inner_radius = center
            outer_radius = top
            orientation = "right"
        elif strand == -1:
            inner_radius = bottom
            outer_radius = center
            orientation = "right"
        else:
            inner_radius = bottom
            outer_radius = top
            orientation = "right"  # backwards compatibility
        return self._draw_arc_arrow(
            inner_radius,
            outer_radius,
            startangle,
            endangle,
            orientation=orientation,
            **kwargs
        )

record = SeqIO.read("Genome.gb", "genbank")
# Setup Data
gd_diagram = GenomeDiagram.Diagram(record.id)
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_track_for_name = gd_diagram.new_track(4, name="name")
gd_feature_set = gd_track_for_features.new_set()
gd_name_set = gd_track_for_name.new_set()

cnt = 0
rate = [0, 0, 0]

#Add Feature to Outcircle

for feature in record.features:
    if feature.type != "gene":
        continue
    
    rate[cnt % 3] += 0.7
    cnt += 1

    color = colors.Color(rate[0], rate[1], rate[2], 1)
    
    gd_name_set.add_feature(
        feature, sigil="ARROW", color=color, label=True, label_size=14, label_angle=0, label_position = "end", arrowshaft_height=0.4
    )

#add Feature to incircle

for feature in record.features:
    if feature.type != "gene":
        continue
    
    color = colors.black

    temp_end = str(feature.location.end)
    temp_start = str(feature.location.start)

    if feature.location.strand < 0:
        temp = temp_end
        temp_end = temp_start
        temp_start = temp

    feature.qualifiers["gene"][0] = temp_end
    gd_feature_set.add_feature(
        feature, sigil="BIGARROW", color=color, label=True, label_size=10, label_angle=0, label_position = "end", arrowshaft_height=0.001,  arrowhead_length=0.0001
    )

    feature.qualifiers["gene"][0] = temp_start
    gd_feature_set.add_feature(
        feature, sigil="BIGARROW", color=color, label=True, label_size=10, label_angle=0, label_position = "start", arrowshaft_height=0.001,  arrowhead_length=0.0001
    )
#Drawing
gd_diagram.draw(
    format="circular",
    circular=True,
    pagesize=(20 * cm, 20 * cm),
    start=0,
    end=len(record),
    circle_core=0.5,
)
gd_diagram.write("circular_genome_graph.png", "PNG")