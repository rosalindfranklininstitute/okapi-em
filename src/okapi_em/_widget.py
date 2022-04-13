"""
This module is an example of a barebones QWidget plugin for napari

It implements the Widget specification.
see: https://napari.org/plugins/stable/guides.html#widgets

Replace code below according to your needs.
"""
#from qtpy.QtWidgets import QWidget, QHBoxLayout, QPushButton, QSlider, QSpinBox, QVBoxLayout, QLabel
#from qtpy import QWidgets
#from qtpy.QtWidgets import *
from turtle import update
from qtpy import QtWidgets
from qtpy.QtCore import Qt
#from napari._qt.widgets.qt_range_slider import QHRangeSlider
#from superqt import QRangeSlider

from . import filters
from . import measure_charging

class MainQWidget(QtWidgets.QWidget):
    # Sets the main window of the widget
    # laying out all controls and doing relevant connections between elements and code

    curData = None
    curDataSlice = None
    myImageLayer = None


    def __init__(self, napari_viewer: 'napari.viewer.Viewer'):
        super().__init__()
        self.viewer = napari_viewer #Reference to the viewer that will be needed later

        #layout = QVBoxLayout()
        self.setLayout(QtWidgets.QVBoxLayout())

        #Put a tab widget that will help select operations
        self.tabwidget = QtWidgets.QTabWidget()
        self.layout().addWidget(self.tabwidget)

        #Create 2 tabs
        self.tab0 = QtWidgets.QWidget()
        self.tabwidget.addTab(self.tab0, "C-Suppression")

        self.tab1 = QtWidgets.QWidget()
        self.tabwidget.addTab(self.tab1, "Tab whatever")

        # /Tab0 \
        self.t0_vbox0 = QtWidgets.QVBoxLayout() #Items organised vertically
        # self.tab0.addWidget( self.t0_vbox0 )
        self.tab0.setLayout(self.t0_vbox0 )

        # Data source --|
        #Group here the controls to select data/slice
        self.group_box0 = QtWidgets.QGroupBox("Data source")
        self.t0_vbox0.addWidget(self.group_box0)

        self.group_box0_vbox0 = QtWidgets.QVBoxLayout()
        self.group_box0.setLayout(self.group_box0_vbox0)

        self.label2 = QtWidgets.QLabel("Set data source")
        self.group_box0_vbox0.addWidget(self.label2)

        self.btnSetDataSourceCurrent = QtWidgets.QPushButton("Set current")
        self.group_box0_vbox0.addWidget(self.btnSetDataSourceCurrent)
        self.btnSetDataSourceCurrent.clicked.connect(self.btnSetDataSourceCurrent_on_click) #Signal
        
        #Grid of labels to show information of data and slice selected
        self.t0_grid0 = QtWidgets.QGridLayout()
        self.group_box0_vbox0.addLayout( self.t0_grid0)

        self.label3 = QtWidgets.QLabel("layer:")
        self.t0_grid0.addWidget(self.label3, 0,0)
        self.labelSelLayerName = QtWidgets.QLabel("")
        self.t0_grid0.addWidget(self.labelSelLayerName, 0,1)
        self.label4 = QtWidgets.QLabel("nslice:")
        self.t0_grid0.addWidget(self.label4, 1,0)
        self.labelSelNSlice = QtWidgets.QLabel("")
        self.t0_grid0.addWidget(self.labelSelNSlice, 1,1)
        self.label5 = QtWidgets.QLabel("slicing axis:")
        self.t0_grid0.addWidget(self.label5, 2,0)
        self.labelSlicingAxis = QtWidgets.QLabel("")
        self.t0_grid0.addWidget(self.labelSlicingAxis, 2,1)
        self.label6 = QtWidgets.QLabel("data slice shape:")
        self.t0_grid0.addWidget(self.label6, 3,0)
        self.labelSliceDims = QtWidgets.QLabel("")
        self.t0_grid0.addWidget(self.labelSliceDims, 3,1)

        #TODO make this a checkbox followed by a box with controls for this filter
        self.label0 = QtWidgets.QLabel("Directional band-pass FFT filter")
        self.t0_vbox0.addWidget( self.label0)

        self.label1 = QtWidgets.QLabel("Select cutoff frequencies")
        self.t0_vbox0.addWidget(self.label1)

        self.t0_hbox = QtWidgets.QHBoxLayout() #Horizontal layout to set low and high freq of filters
        self.t0_vbox0.addLayout( self.t0_hbox)
        
        self.label_low = QtWidgets.QLabel("low:")
        self.t0_hbox.addWidget(self.label_low)

        self.freqselector_low = QtWidgets.QSpinBox()
        self.freqselector_low.setMinimum(0)
        self.freqselector_low.setMaximum(100)
        self.freqselector_low.setValue(0)
        self.freqselector_low.setSingleStep(1)
        self.t0_hbox.addWidget(self.freqselector_low)
        self.freqselector_low.valueChanged.connect(self.freqselector_low_valueChanged)

        self.label_high = QtWidgets.QLabel("high:")
        self.t0_hbox.addWidget(self.label_high)

        self.freqselector_high = QtWidgets.QSpinBox()
        self.freqselector_high.setMinimum(0)
        self.freqselector_high.setMaximum(100)
        self.freqselector_high.setValue(5)
        self.freqselector_high.setSingleStep(1)
        print(f"self.freqselector_high.stepType: {self.freqselector_high.stepType()}")
        #self.freqselector_high.setStepType(QtWidgets.QAbstractSpinBox.)
        self.t0_hbox.addWidget(self.freqselector_high)
        self.freqselector_high.valueChanged.connect(self.freqselector_high_valueChanged)

        # self.layout().addWidget(self.label1)
        # self.layout().addWidget(self.freqselector)

        # self.freqRangeSlider = QRangeSlider(Qt.Horizontal)
        # self.freqRangeSlider.setRange(0,100)
        # self.freqRangeSlider.setValue((0,5))
        # self.t0_vbox0.addWidget(self.freqRangeSlider)

        self.chkPreview = QtWidgets.QCheckBox("Preview")
        self.chkPreview.setChecked(True)
        self.t0_vbox0.addWidget(self.chkPreview)
        #self.chkPreview.stateChanged.connect(lambda:self.updatePreview(self.chkPreview)) #Signal
        self.chkPreview.toggled.connect(self.updatePreview)

        self.btnApply = QtWidgets.QPushButton("Apply")
        self.t0_vbox0.addWidget(self.btnApply)
        self.btnApply.clicked.connect(self.btnApply_on_click) #Signal


    def freqselector_high_valueChanged(self):
        #check value is not lower than lowvalue
        v0 = self.freqselector_high.value()
        if v0< self.freqselector_low.value():
            self.freqselector_high.setValue( v0+1)
        self.updatePreview()
    
    def freqselector_low_valueChanged(self):
        v0 = self.freqselector_low.value()
        if v0> self.freqselector_high.value():
            self.freqselector_low.setValue(v0-1)
        self.updatePreview()

    def btnApply_on_click(self):
        #From current layer/data selected, apply the filter to the whole data
        print("NOT IMPLEMENTED YET")
        #TODO

    def updatePreview(self):
        #Check if the preview check box is enabled
        if self.chkPreview.isChecked():
            #Check if there is valid data for creating the preview
            
            if not self.curDataSlice is None:
                #Debug
                print(f"curData.shape = {self.curData.shape}")
                self.labelSliceDims.setText(str(self.curDataSlice.shape))

                freq_low = self.freqselector_low.value()
                freq_high = self.freqselector_high.value()

                print(f"freq_low : {freq_low} , freq_high: {freq_high}")

                #Calculate the filtered data
                dataPreview = filters.fft_bandpass_filter_dirx_2D(self.curDataSlice,freq_low, freq_high)

                if self.myImageLayer is None or not self.myImageLayer in self.viewer.layers:
                    self.myImageLayer = self.viewer.add_image(dataPreview, name="preview")
                    self.myImageLayer._keep_auto_contrast=True #This helps show image with auto contrast however associated button does not show it is pressed
                    #self.myImageLayer.reset_contrast_limits() #Try, doesnt work
                else:
                    #Dont create a new layer, just modify it
                    self.myImageLayer.data = dataPreview 
                    #self.myImageLayer.refresh()

    def btnSetDataSourceCurrent_on_click(self):
        #Check if there is valid data for creating the preview
        
        #Get current layer and current slice if applicable
        active0 = self.viewer.layers.selection.active
        if not active0 is None:
            print("Current active layer name is: ", active0.name)
            self.labelSelLayerName.setText(active0.name)

            self.curData = active0.data

            self.curDataSlice  = None #For the previews
            if self.curData.ndim==3:
                print("ndim=3")

                #Get the current slice if the viewer is in 2D mode

                if self.viewer.dims.ndisplay == 2:
                    #Get the currently viewing slice
                    dims0 = self.viewer.dims

                    nslice = dims0.current_step[0] #Gets the slice number that is currently being viewed, 1st index is the slice number, but which axis is not given here
                    print("nslice: ", str(nslice))
                    self.labelSelNSlice.setText(str(nslice))

                    order0 = dims0.order[0] #Gets the axis (or plane), that the Viewer is currently viewing
                    if order0 == 0: #z
                        self.curDataSlice = self.curData[nslice,...]
                        self.labelSlicingAxis.setText("z")
                    elif order0 == 1: #y
                        self.curDataSlice = self.curData[:,nslice,:]
                        self.labelSlicingAxis.setText("y")
                    else: #x
                        self.curDataSlice = self.curData[...,nslice]
                        self.labelSlicingAxis.setText("x")

                    #Get the max freq for the filter from the curDataSlice
                    maxfreq = int(self.curDataSlice.shape[1])
                    #Adjust max and min values of freq range numerics
                    self.freqselector_low.setMaximum(maxfreq)
                    self.freqselector_high.setMaximum(maxfreq)

                    self.updatePreview()

                else:
                    print("Viewer is in 3D mode, cannot get slice")
                    self.labelSelNSlice.setText("Viewer in 3D mode, cannot get slice")

            elif self.curData.ndim==2:
                print("ndim=2")
                self.curDataSlice = self.curData
                self.updatePreview()
            else: #ndims=1
                print("ndim=1, Unsupported")
                return
        else:
            print("No active layer to select data source.")
