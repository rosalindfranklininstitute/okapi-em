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
        self.tabwidget.addTab(self.tab1, "Measure")

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

        # Radiobutton followed by a box with controls for this filter
        self.rbDirFFT = QtWidgets.QRadioButton("Directional band-pass FFT filter")
        self.rbDirFFT.setChecked(True)
        self.t0_vbox0.addWidget( self.rbDirFFT)
        #TODO, may need 'connection'

        self.grpBox1 = QtWidgets.QGroupBox("") #Box, no title, radiobutton does it
        self.t0_vbox0.addWidget(self.grpBox1)
        self.vbox2 = QtWidgets.QVBoxLayout()
        self.grpBox1.setLayout(self.vbox2)

        self.label1 = QtWidgets.QLabel("Select cutoff frequencies")
        self.vbox2.addWidget(self.label1)

        self.t0_hbox = QtWidgets.QHBoxLayout() #Horizontal layout to set low and high freq of filters
        self.vbox2.addLayout( self.t0_hbox)
        
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
        self.t0_hbox.addWidget(self.freqselector_high)
        self.freqselector_high.valueChanged.connect(self.freqselector_high_valueChanged)

        self.chkPreview = QtWidgets.QCheckBox("Preview")
        self.chkPreview.setChecked(True)
        self.t0_vbox0.addWidget(self.chkPreview)
        #self.chkPreview.stateChanged.connect(lambda:self.updatePreview(self.chkPreview)) #Signal
        self.chkPreview.toggled.connect(self.updatePreview)

        self.btnApply = QtWidgets.QPushButton("Apply")
        self.t0_vbox0.addWidget(self.btnApply)
        self.btnApply.clicked.connect(self.btnApply_on_click) #Signal


        # /Tab1 \
        self.vBox3 = QtWidgets.QVBoxLayout() #Items organised vertically
        self.tab1.setLayout(self.vBox3 )

        #Add stuff
        self.grpBox2 = QtWidgets.QGroupBox("Measure charging")
        self.vBox3.addWidget(self.grpBox2)
        self.grpBox2_vlayout = QtWidgets.QVBoxLayout()
        self.grpBox2.setLayout(self.grpBox2_vlayout)

        self.lbl7 = QtWidgets.QLabel("Please ensure correct data source in layer list is selected before calculation")
        self.lbl7.setWordWrap(True)
        self.grpBox2_vlayout.addWidget(self.lbl7)

        self.hBox2 = QtWidgets.QHBoxLayout()
        self.grpBox2_vlayout.addLayout( self.hBox2)
        self.lblMCTileSize = QtWidgets.QLabel("Tile Size")
        self.hBox2.addWidget(self.lblMCTileSize)
        self.mc_tilesize_selector = QtWidgets.QSpinBox()
        self.mc_tilesize_selector.setMinimum(8)
        self.mc_tilesize_selector.setMaximum(2048)
        self.mc_tilesize_selector.setValue(256)
        self.hBox2.addWidget(self.mc_tilesize_selector)

        self.btnMCCalculate = QtWidgets.QPushButton("Calculate")
        self.grpBox2_vlayout.addWidget(self.btnMCCalculate)
        self.btnMCCalculate.clicked.connect(self.btnMCCalculate_onclick) #Signal




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
    
    def btnMCCalculate_onclick(self):
        #Run the measure charging

        #Read tile size parameter and gets current image
        tilesize = self.mc_tilesize_selector.value()

        active0 = self.viewer.layers.selection.active
        data0= active0.data

        curDataSlice=None
        if data0.ndim==2:
            curDataSlice=data0
        elif data0.ndim==3:
            dims0 = self.viewer.dims
            nslice = dims0.current_step[0] #Gets the slice number that is currently being viewed, 1st index is the slice number, but which axis is not given here
            orderdict = { 0:'z' , 1:'y' , 2:'x'}
            order0 = dims0.order[0]
            print(f"nslice: {str(nslice)} ; order: {orderdict[order0]}")

            if order0 == 0: #z
                curDataSlice = data0[nslice,...]
            elif order0 == 1: #y
                curDataSlice = data0[:,nslice,:]
            else: #x
                curDataSlice = data0[...,nslice]
        
        if not curDataSlice is None:
            #Calculate
            overlay = measure_charging.generate_heatmap(curDataSlice, tilesize)
            self.viewer.add_image(overlay, name="Charging artefacts", opacity=0.3, colormap="viridis")


