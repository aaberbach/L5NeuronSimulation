import streamlit as st
import pandas as pd
from numpy.random import random
import plotly.graph_objects as go

segs = pd.read_csv('L5Segments.csv')
segs_degrees = pd.read_csv('SegmentsDegrees.csv').groupby(['Type','Sec ID'])['Degrees'].max().reset_index()
segs['segmentID'] = segs.index
segs = segs.set_index(['Type','Sec ID']).join(segs_degrees.set_index(['Type','Sec ID'])).reset_index()


trace1 = [go.Scatter3d(x=segs['Coord X'], 
                       y=segs['Coord Z'], 
                       z=segs['Coord Y'], 
                       mode='markers',
                       opacity=0.2)]

fig = go.Figure(trace1, layout=go.Layout(width=10,height=2000))

st.plotly_chart(fig, use_container_width=True)
