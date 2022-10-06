"""
This module is an example of a barebones QWidget plugin for napari

It implements the Widget specification.
see: https://napari.org/plugins/stable/guides.html#widgets

Replace code below according to your needs.
"""
#from qtpy.QtWidgets import QWidget, QHBoxLayout, QPushButton, QSlider, QSpinBox, QVBoxLayout, QLabel
#from qtpy import QWidgets
from qtpy import QtWidgets

#from qtpy.QtCore import Qt
#from napari._qt.widgets.qt_range_slider import QHRangeSlider
#from superqt import QRangeSlider

import napari
from napari.types import ImageData

from . import filters
from . import measure_charging

#QtWidgets.QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)

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

        #self.t0_vbox0 = QtWidgets.QVBoxLayout() #Items organised vertically
        # self.tabChargeSuppr.addWidget( self.t0_vbox0 )
        #self.tabChargeSuppr.setLayout(self.t0_vbox0 )
        #self.layout().addWidget(self.t0_vbox0)

        # Data source Box |---|
        #Group here the controls to select data/slice
        self.group_box0 = QtWidgets.QGroupBox("Data source")
        #self.t0_vbox0.addWidget(self.group_box0)
        self.layout().addWidget(self.group_box0)

        self.group_box0_vbox0 = QtWidgets.QVBoxLayout()
        self.group_box0.setLayout(self.group_box0_vbox0)

        self.label2 = QtWidgets.QLabel("Set data source")
        self.group_box0_vbox0.addWidget(self.label2)
        # self.btnSetDataSourceCurrent = QtWidgets.QPushButton("Set currently selected")
        # self.group_box0_vbox0.addWidget(self.btnSetDataSourceCurrent)
        # self.btnSetDataSourceCurrent.clicked.connect(self.btnSetDataSourceCurrent_on_click) #Signal
        
        self.w_setselect = widget_ImageSetCurrrentSelect(napari_viewer)
        self.group_box0_vbox0.addWidget(self.w_setselect)

        # #Grid of labels to show information of data and slice selected
        # self.t0_grid0 = QtWidgets.QGridLayout()
        # self.group_box0_vbox0.addLayout( self.t0_grid0)

        # self.label3 = QtWidgets.QLabel("layer:")
        # self.t0_grid0.addWidget(self.label3, 0,0)
        # self.labelSelLayerName = QtWidgets.QLabel("")
        # self.t0_grid0.addWidget(self.labelSelLayerName, 0,1)
        # self.label4 = QtWidgets.QLabel("nslice:")
        # self.t0_grid0.addWidget(self.label4, 1,0)
        # self.labelSelNSlice = QtWidgets.QLabel("")
        # self.t0_grid0.addWidget(self.labelSelNSlice, 1,1)
        # self.label5 = QtWidgets.QLabel("slicing axis:")
        # self.t0_grid0.addWidget(self.label5, 2,0)
        # self.labelSlicingAxis = QtWidgets.QLabel("")
        # self.t0_grid0.addWidget(self.labelSlicingAxis, 2,1)
        # self.label6 = QtWidgets.QLabel("data slice shape:")
        # self.t0_grid0.addWidget(self.label6, 3,0)
        # self.labelSliceDims = QtWidgets.QLabel("")
        # self.t0_grid0.addWidget(self.labelSliceDims, 3,1)
        # Data Source Box |___| end

        #Put a tab widget that will help select operations
        self.tabwidget = QtWidgets.QTabWidget()
        self.layout().addWidget(self.tabwidget)
        # and set its layout as vertical
        #self.t0_vbox0.addWidget(self.tabwidget)

        #Create all tabs here
        
        # /tabChargeSuppr \ t0_vbox0
        self.tabChargeSuppr = QtWidgets.QWidget()
        self.tabwidget.addTab(self.tabChargeSuppr, "C-Suppression")
        # Tab with charge suppression filters
        # Contains FFT directional filter plus others

        # /tabMeasChargingArtif \
        self.tabMeasChargingArtif = QtWidgets.QWidget()
        self.tabwidget.addTab(self.tabMeasChargingArtif, "Measure") #Measure charging artifacts

        # /tabQuoll \
        self.tabQuoll = QtWidgets.QWidget()
        self.tabwidget.addTab(self.tabQuoll, "Quoll")


        self.t0_vbox0 = QtWidgets.QVBoxLayout() #Items organised vertically
        self.tabChargeSuppr.setLayout(self.t0_vbox0)

        # Filter: Directional band-pass
        # Radiobutton followed by a box with controls for this filter
        self.rbDirFFT = QtWidgets.QRadioButton("Directional band-pass FFT filter")
        self.rbDirFFT.setChecked(True)
        self.t0_vbox0.addWidget( self.rbDirFFT)
        #TODO, may need 'connection'
        #Options for this filter are contained in the next box grpBox1
        self.grpBox1 = QtWidgets.QGroupBox("") #Box, no title, radiobutton does it
        #self.tabChargeSuppr.addWidget(self.grpBox1)
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

        self.rbChaffer = QtWidgets.QRadioButton("Chaffer")
        self.rbChaffer.setChecked(False)
        self.t0_vbox0.addWidget(self.rbChaffer)
        self.grpBoxChaffer = QtWidgets.QGroupBox("")
        self.t0_vbox0.addWidget(self.grpBoxChaffer) #Add chaffer UI to this grpBoxChaffer
        #Need select source file,
        #Need user to select labels file


        self.btnFFTDirFilterSlice = QtWidgets.QPushButton("Filter (slice)")
        self.t0_vbox0.addWidget(self.btnFFTDirFilterSlice)
        self.btnFFTDirFilterSlice.clicked.connect(self.doFFTDirFilterSlice)

        self.btnApply = QtWidgets.QPushButton("Filter")
        self.t0_vbox0.addWidget(self.btnApply)
        self.btnApply.clicked.connect(self.btnApply_on_click) #Signal



        self.vBox3 = QtWidgets.QVBoxLayout() #Items organised vertically
        self.tabMeasChargingArtif.setLayout(self.vBox3 )

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
        self.doFFTDirFilterSlice()
    
    def freqselector_low_valueChanged(self):
        v0 = self.freqselector_low.value()
        if v0> self.freqselector_high.value():
            self.freqselector_low.setValue(v0-1)
        self.doFFTDirFilterSlice()

    def btnApply_on_click(self):
        #From current layer/data selected, apply the filter to the whole data
        print("NOT IMPLEMENTED YET")
        #TODO

    def doFFTDirFilterSlice(self):
        #Check if the preview check box is enabled
        if self.chkPreview.isChecked():
            #Check if there is valid data for creating the preview
            
            if not self.curDataSlice is None:
                #Debug
                #print(f"curData.shape = {self.curData.shape}")
                self.labelSliceDims.setText(str(self.curDataSlice.shape))

                freq_low = self.freqselector_low.value()
                freq_high = self.freqselector_high.value()

                #print(f"freq_low : {freq_low} , freq_high: {freq_high}")

                #Calculate the filtered data
                dataPreview = filters.fft_bandpass_filter_dirx_2D(self.curDataSlice,freq_low, freq_high)

                if self.myImageLayer is None or not self.myImageLayer in self.viewer.layers:
                    #self.myImageLayer = self.viewer.add_image(dataPreview, name="preview")
                    #self.myImageLayer._keep_auto_contrast=True #This helps show image with auto contrast however associated button does not show it is pressed
                    #self.myImageLayer.reset_contrast_limits() #Try, doesnt work
                    self.myImageLayer = self.viewer.add_image(dataPreview, name="preview", contrast_limits= (self.curDataSlice.min(), self.curDataSlice.max()))

                else:
                    #Dont create a new layer, just modify it
                    self.myImageLayer.data = dataPreview 
                    #self.myImageLayer.refresh()

    # def btnSetDataSourceCurrent_on_click(self):
    #     #Check if there is valid data for creating the preview
        
    #     #Get current layer and current slice if applicable
    #     active0 = self.viewer.layers.selection.active
    #     if not active0 is None:
    #         print("Current active layer name is: ", active0.name)
    #         self.labelSelLayerName.setText(active0.name)

    #         self.curData = active0.data

    #         self.curDataSlice  = None #For the previews
    #         if self.curData.ndim==3:
    #             print("ndim=3")

    #             #Get the current slice if the viewer is in 2D mode

    #             if self.viewer.dims.ndisplay == 2:
    #                 #Get the currently viewing slice
    #                 dims0 = self.viewer.dims

    #                 nslice = dims0.current_step[0] #Gets the slice number that is currently being viewed, 1st index is the slice number, but which axis is not given here
    #                 print("nslice: ", str(nslice))
    #                 self.labelSelNSlice.setText(str(nslice))

    #                 order0 = dims0.order[0] #Gets the axis (or plane), that the Viewer is currently viewing
    #                 if order0 == 0: #z
    #                     self.curDataSlice = self.curData[nslice,...]
    #                     self.labelSlicingAxis.setText("z")
    #                 elif order0 == 1: #y
    #                     self.curDataSlice = self.curData[:,nslice,:]
    #                     self.labelSlicingAxis.setText("y")
    #                 else: #x
    #                     self.curDataSlice = self.curData[...,nslice]
    #                     self.labelSlicingAxis.setText("x")

    #                 #Get the max freq for the filter from the curDataSlice
    #                 maxfreq = int(self.curDataSlice.shape[1])
    #                 #Adjust max and min values of freq range numerics
    #                 self.freqselector_low.setMaximum(maxfreq)
    #                 self.freqselector_high.setMaximum(maxfreq)

    #                 self.updatePreview()

    #             else:
    #                 print("Viewer is in 3D mode, cannot get slice")
    #                 self.labelSelNSlice.setText("Viewer in 3D mode, cannot get slice")

    #         elif self.curData.ndim==2:
    #             print("ndim=2")
    #             self.curDataSlice = self.curData
    #             self.updatePreview()
    #         else: #ndims=1
    #             print("ndim=1, Unsupported")
    #             return
    #     else:
    #         print("No active layer to select data source.")
    
    def btnMCCalculate_onclick(self):
        #Run the measure charging

        print("btnMCCalculate_onclick()")

        #Read tile size parameter and gets current image
        tilesize = self.mc_tilesize_selector.value()

        data0 = self.w_setselect.get_active_image_data_slice()

        curDataSlice=data0
        if not curDataSlice is None:
            #Calculate
            overlay = measure_charging.generate_heatmap(curDataSlice, tilesize)
            self.viewer.add_image(overlay, name="Charging artefacts", opacity=0.3, colormap="viridis")

    # def get_napari_image_list(self):
    #     #Code similar to found in
    #     # https://github.com/haesleinhuepf/napari-accelerated-pixel-and-object-classification/blob/main/napari_accelerated_pixel_and_object_classification/_dock_widget.py#L392
    #     # for l in self.viewer.layers:
    #     #     if isinstance(l, napari.layers.Image):
    #     #         suffix = ""
    #     #         if len(l.data.shape) == 4:
    #     #             suffix = " (current timepoint)"
    #     #         item = QListWidgetItem(l.name + suffix)
    #     #         self._available_images.append(l)
    #     #         self.image_list.addItem(item)
    #     #         if l in selected_images:
    #     #             item.setSelected(True)
    #     pass



class widget_ImageSetCurrrentSelect(QtWidgets.QWidget):
    '''
    A QtWidget that helps select Image data for further processing

    Needs to be initialised with a napari Viewer object

    '''

    def __init__(self, napari_viewer: 'napari.viewer.Viewer'):
        super().__init__()
        self.viewer = napari_viewer #Reference to the viewer that will be needed later

        layout = QtWidgets.QGridLayout()
        self.setLayout(layout)

        self.b_set = QtWidgets.QPushButton("Set currently selected 2")
        layout.addWidget(self.b_set,0,0)
        self.b_set.clicked.connect(self.b_set_on_click) #Signal
  
        layout.addWidget(QtWidgets.QLabel("name"),1,0)
        self.l_name_v = QtWidgets.QLabel("")
        layout.addWidget(self.l_name_v,1,1)

        layout.addWidget(QtWidgets.QLabel("shape"),2,0)
        self.l_shape_v = QtWidgets.QLabel("")
        layout.addWidget(self.l_shape_v,2,1)

        layout.addWidget(QtWidgets.QLabel("slicing axis"),3,0)
        self.l_slicingaxis_v = QtWidgets.QLabel("")
        layout.addWidget(self.l_slicingaxis_v,3,1)

        layout.addWidget(QtWidgets.QLabel("slice shape"),4,0)
        self.l_sliceshape_v = QtWidgets.QLabel("")
        layout.addWidget(self.l_sliceshape_v,4,1)

        layout.addWidget(QtWidgets.QLabel("nslice"),5,0)
        self.l_nslice_v = QtWidgets.QLabel("")
        layout.addWidget(self.l_nslice_v,5,1)

        self.curImage=None

    def b_set_on_click(self):
        print("b_set_on_click()")
        active0=self.get_active_image() #This also sets default image
        if not active0 is None:
            self.l_name_v.setText(active0.name)

            #Check layer is image
            #self.set_default_napari_image(active0)

    def set_default_napari_image(self,nap_image):
        print("set_default_napari_image()")
        self.curImage=None
        if isinstance(nap_image, napari.layers.Image):
            self.l_name_v.setText(nap_image.name)
            print("Current active layer name is: ", nap_image.name)
            self.curImage = nap_image
            self.l_shape_v.setText(str(self.curImage.data.shape))
            

    def get_active_image(self):
        active0=None
        if not self.viewer is None:
            active0 = self.viewer.layers.selection.active #check if any layer is active
            self.set_default_napari_image(active0)

        return active0

    def get_image_selected_data(self) ->ImageData:
        print("get_image_selected_data()")
        if self.curImage is None:
            self.b_set_on_click() #Does the same thing as clicking the button to set default image data
            return self.curImage.data
        
        return self.curImage.data
    

    def get_active_image_data_slice(self):
        print("get_active_image_data_slice()")

        curDataSlice=None

        curData = self.get_image_selected_data()
        #print(f"curData:{curData}")

        if not curData is None:
            if curData.ndim==3:
                if self.viewer.dims.ndisplay == 2:
                    #Get the currently viewing slice
                    dims0 = self.viewer.dims

                    nslice = dims0.current_step[0] #Gets the slice number that is currently being viewed, 1st index is the slice number, but which axis is not given here
                    print("nslice: ", str(nslice))
                    self.l_nslice_v.setText(str(nslice))

                    order0 = dims0.order[0] #Gets the axis (or plane), that the Viewer is currently viewing
                    if order0 == 0: #z
                        curDataSlice = self.curImage.data[nslice,...]
                        self.l_slicingaxis_v.setText('z')
                    elif order0 == 1: #y
                        curDataSlice = self.curImage.data[:,nslice,:]
                        self.l_slicingaxis_v.setText("y")
                    else: #x
                        curDataSlice = self.curImage.data[...,nslice]
                        self.l_slicingaxis_v.setText("x")

                    #Slice shape
                    s_shape = curDataSlice.shape
                    self.l_sliceshape_v.setText(str(s_shape))

                else:
                    print("Viewer is in 3D mode, cannot get slice")
                    self.l_slicingaxis_v.setText("Viewer in 3D mode, cannot get slice")
            else:
                print("ndim=2. The image layer itself is the slice.")
                curDataSlice= curData
            
        return curDataSlice

