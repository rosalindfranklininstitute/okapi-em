'''
Copyright 2022 Rosalind Franklin Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

'''
Main napari widget with data selector and tabs with the tools
'''

#from qtpy.QtWidgets import QWidget, QHBoxLayout, QPushButton, QSlider, QSpinBox, QVBoxLayout, QLabel
#from qtpy import QWidgets
from qtpy import QtWidgets

#from qtpy.QtCore import Qt
#from napari._qt.widgets.qt_range_slider import QHRangeSlider
#from superqt import QRangeSlider

import napari
from napari.types import ImageData, LayerDataTuple
from napari.utils import progress

from . import filters
from . import measure_charging

import numpy as np

# from magicgui import magicgui
# from magicgui import magic_factory

import chafer
from . import quoll_wrapper as quoll


from . import slice_alignment

class MainQWidget(QtWidgets.QWidget):
    # Sets the main window of the widget
    # laying out all controls and doing relevant connections between elements and code

    curData = None
    curDataSlice = None
    myFFTFilterDirSliceRes = None


    def __init__(self, napari_viewer: 'napari.viewer.Viewer'):
        super().__init__()
        self.viewer = napari_viewer #Reference to the viewer that will be needed later

        print(f"napari_viewer: {napari_viewer}")

        #layout = QVBoxLayout()
        self.setLayout(QtWidgets.QVBoxLayout())

        #self.t0_vbox0 = QtWidgets.QVBoxLayout() #Items organised vertically
        # self.tabChargeSuppr.addWidget( self.t0_vbox0 )
        #self.tabChargeSuppr.setLayout(self.t0_vbox0 )
        #self.layout().addWidget(self.t0_vbox0)

        # Data source Box |---|
        #Group here the controls to select data/slice
        # self.group_box0 = QtWidgets.QGroupBox("Data source")
        # #self.t0_vbox0.addWidget(self.group_box0)
        # self.layout().addWidget(self.group_box0)

        # self.group_box0_vbox0 = QtWidgets.QVBoxLayout()
        # self.group_box0.setLayout(self.group_box0_vbox0)

        # self.group_box0_vbox0.addWidget(QtWidgets.QLabel("Set data source"))
        # # self.btnSetDataSourceCurrent = QtWidgets.QPushButton("Set currently selected")
        # # self.group_box0_vbox0.addWidget(self.btnSetDataSourceCurrent)
        # # self.btnSetDataSourceCurrent.clicked.connect(self.btnSetDataSourceCurrent_on_click) #Signal
        
        # self.w_setselect = widget_ImageSetCurrrentSelect(napari_viewer)
        # self.group_box0_vbox0.addWidget(self.w_setselect)
        self.w_setselect = widget_ImageSetCurrentSelect(napari_viewer, "Data source")
        self.layout().addWidget(self.w_setselect)

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
        self.ui_do_tab_ChargeSuppression()
        self.ui_do_tab_SliceAlignment()
        self.ui_do_tab_Quoll()

    def ui_do_tab_ChargeSuppression(self):
        # /tabChargeSuppr \ t0_vbox0
        self.tabChargeSuppr = QtWidgets.QWidget()
        self.tabwidget.addTab(self.tabChargeSuppr, "C-Suppression")
        # Tab with charge suppression filters
        # Contains FFT directional filter plus others

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
        self.grpBoxChaffer_layout = QtWidgets.QVBoxLayout()
        self.grpBoxChaffer.setLayout(self.grpBoxChaffer_layout)

        #label data selector
        self.w_setselectlabel = widget_ImageSetCurrentSelect(napari_viewer, "Labels data source")
        self.grpBoxChaffer_layout.addWidget(self.w_setselectlabel)

        #need parameters
        #nlinesaverage=20,
        #data_fit_max_length_edge_px = 700,
        #data_fit_min_length_px = 50

        qgrid2_layout = QtWidgets.QGridLayout()
        self.grpBoxChaffer_layout.addLayout(qgrid2_layout)

        qgrid2_layout.addWidget( QtWidgets.QLabel("nlinesaverage") ,0,0)
        self.chafer_widg_nlines = QtWidgets.QSpinBox()
        self.chafer_widg_nlines.setMinimum(1)
        self.chafer_widg_nlines.setMaximum(65536)
        self.chafer_widg_nlines.setValue(20)
        qgrid2_layout.addWidget( self.chafer_widg_nlines ,0,1)

        qgrid2_layout.addWidget( QtWidgets.QLabel("fit length min px") ,1,0)
        self.chafer_widg_minlength = QtWidgets.QSpinBox()
        self.chafer_widg_minlength.setMinimum(1)
        self.chafer_widg_minlength.setMaximum(65536)
        self.chafer_widg_minlength.setValue(50)
        qgrid2_layout.addWidget( self.chafer_widg_minlength ,1,1)

        qgrid2_layout.addWidget( QtWidgets.QLabel("fit length max px") ,2,0)
        self.chafer_widg_maxlength = QtWidgets.QSpinBox()
        self.chafer_widg_maxlength.setMinimum(1)
        self.chafer_widg_maxlength.setMaximum(65536)
        self.chafer_widg_maxlength.setValue(700)
        qgrid2_layout.addWidget( self.chafer_widg_maxlength ,2,1)
        
        #not working
        # self.chafer_widget = chaffer_magicgui_widget.native
        # self.t0_vbox0.addWidget(self.chafer_widget)
        #self.t0_vbox0.addWidget(chafer_widget.native)

        #Need select source file,
        #Need user to select labels file

        self.btnFFTDirFilterSlice = QtWidgets.QPushButton("Filter (slice)")
        self.t0_vbox0.addWidget(self.btnFFTDirFilterSlice)
        self.btnFFTDirFilterSlice.clicked.connect(self.doFFTDirFilterSlice)

        self.btnApply = QtWidgets.QPushButton("Filter")
        self.t0_vbox0.addWidget(self.btnApply)
        self.btnApply.clicked.connect(self.doFFTDirFilterWhole) #Signal

    def ui_do_tab_SliceAlignment(self):
        # /tabSliceAlignemnt \
        self.tabSliceAlignment = QtWidgets.QWidget()
        self.tabwidget.addTab(self.tabSliceAlignment, "Slice Alignment")
        self.tabSliceAlignment_l = QtWidgets.QVBoxLayout() #Items organised vertically
        self.tabSliceAlignment.setLayout(self.tabSliceAlignment_l)
        #Button to start the alignement
        self.btnSACalculate = QtWidgets.QPushButton("Align")
        self.tabSliceAlignment_l.addWidget(self.btnSACalculate)
        self.btnSACalculate.clicked.connect(self.btnSACalculate_onclick) #Signal

    def ui_do_tab_Quoll(self):
        # /tabQuoll \
        if quoll.is_quoll_available:
            self.tabQuoll = QtWidgets.QWidget()
            self.tabwidget.addTab(self.tabQuoll, "Quoll")

            self.vBox4 = QtWidgets.QVBoxLayout() #Items organised vertically
            self.tabQuoll.setLayout(self.vBox4 )

            # self.vBox4.addWidget(QtWidgets.QLabel("FRC calculations"))
            # self.btnQuollCalcFRC = QtWidgets.QPushButton("Calculate FRC (whole slice)")
            # self.vBox4.addWidget(self.btnQuollCalcFRC)
            # self.btnQuollCalcFRC.clicked.connect(self.btnQuollCalcFRC_onclick) #Signal

            # self.hBox4 = QtWidgets.QHBoxLayout()
            # self.vBox4.addLayout(self.hBox4)
            # self.hBox4.addWidget(QtWidgets.QLabel("Value /px"))
            # self.lblFRCWholeImg=QtWidgets.QLabel("0")
            # self.hBox4.addWidget(self.lblFRCWholeImg)

            self.btnQuollCalcFRCTiled = QtWidgets.QPushButton("Calculate FRC tiled")
            self.vBox4.addWidget(self.btnQuollCalcFRCTiled)
            self.btnQuollCalcFRCTiled.clicked.connect(self.btnQuollCalcFRCTiled_onclick) #Signal
            qgrid3_layout = QtWidgets.QGridLayout()
            self.vBox4.addLayout(qgrid3_layout)
            qgrid3_layout.addWidget(QtWidgets.QLabel("tile size (px):"),0,0)
            self.spbQuolTileSize = QtWidgets.QSpinBox()
            self.spbQuolTileSize.setMinimum(8)
            self.spbQuolTileSize.setMaximum(10000)
            self.spbQuolTileSize.setValue(256)
            qgrid3_layout.addWidget(self.spbQuolTileSize, 0,1)
            self.quollTableContainerWidget=QtWidgets.QWidget()
            qgrid3_layout.addWidget(self.quollTableContainerWidget)
            
            #table with infomation
            self.quollTableContainerWidgetLayout=QtWidgets.QVBoxLayout()
            self.quollTableContainerWidget.setLayout(self.quollTableContainerWidgetLayout)
            self.quollTable=widget_PandasDFTable(self.viewer)
            self.quollTableContainerWidgetLayout.addWidget(self.quollTable)


    def freqselector_high_valueChanged(self):
        #check value is not lower than lowvalue
        v0 = self.freqselector_high.value()
        if v0< self.freqselector_low.value():
            self.freqselector_high.setValue( 0)
        #self.doFFTDirFilterSlice()
    
    def freqselector_low_valueChanged(self):
        v0 = self.freqselector_low.value()
        if v0> self.freqselector_high.value():
            self.freqselector_low.setValue(v0-1)
        #self.doFFTDirFilterSlice()

    def doFFTDirFilterSlice(self):
        
        print("doFFTDirFilterSlice()")

        curDataSlice = self.w_setselect.get_active_image_selected_data_slice()
        if not curDataSlice is None:

            # myfunc=None

            # freq_low = self.freqselector_low.value()
            # freq_high = self.freqselector_high.value()
            # def func1(data2D):
            #     return filters.fft_bandpass_filter_dirx_2D(data2D,freq_low, freq_high)

            # c_filter = chafer.cls_charge_artifact_suppression_filter(
            #     self.chafer_widg_nlines.value(),
            #     self.chafer_widg_maxlength.value(),
            #     self.chafer_widg_minlength.value())

            # def func2(data2d, dlabels):
            #     res,opts= c_filter.charge_artifact_FD_filter_downup_av_prevlines3_2d(data2d,dlabels)
            #     return res

            datares=None

            #Check which filter to use
            if self.rbDirFFT.isChecked():
                freq_low = self.freqselector_low.value()
                freq_high = self.freqselector_high.value()
                
                datares = filters.fft_bandpass_filter_dirx_2D(curDataSlice,freq_low, freq_high)

            elif self.rbChaffer.isChecked():
                curLabelsSlice= self.w_setselectlabel.get_active_image_selected_data_slice()
                c_filter = chafer.cls_charge_artifact_suppression_filter(
                    self.chafer_widg_nlines.value(),
                    self.chafer_widg_maxlength.value(),
                    self.chafer_widg_minlength.value()
                    )
                
                datares, opts = c_filter.charge_artifact_FD_filter_downup_av_prevlines3_2d(curDataSlice, curLabelsSlice)

            if not datares is None:
                #Not sure if using a single preview data layer is ok
                if self.myFFTFilterDirSliceRes is None or not self.myFFTFilterDirSliceRes in self.viewer.layers:
                    self.myFFTFilterDirSliceRes = self.viewer.add_image(datares, name="c-filtered", contrast_limits= (curDataSlice.min(), curDataSlice.max()))
                else:
                    #Dont create a new layer, just modify it
                    self.myFFTFilterDirSliceRes.data = datares
                    #self.myImageLayer.refresh()

    def doFFTDirFilterWhole(self):
        #From current layer/data selected, apply the filter to the whole data
        
        print("doFFTDirFilterWhole()")
        curData = self.w_setselect.get_active_image_selected_data()

        if not curData is None:
            
            data_filt_whole=None
            
            cCalculate=None

            #Check which filter to use
            if self.rbDirFFT.isChecked():
                freq_low = self.freqselector_low.value()
                freq_high = self.freqselector_high.value()
            
                def func(data2d):
                    return filters.fft_bandpass_filter_dirx_2D(data2d,freq_low, freq_high)

                cCalculate = cSliceBySliceFunctHelper(self.viewer,func, curData)

            
            elif self.rbChaffer.isChecked():
                #Dont use the chaffer's 3d method
                #Use this engine to do slice by slice so that it can be monitored

                curLabels= self.w_setselectlabel.get_active_image_selected_data()

                c_filter = chafer.cls_charge_artifact_suppression_filter(
                    self.chafer_widg_nlines.value(),
                    self.chafer_widg_maxlength.value(),
                    self.chafer_widg_minlength.value()
                    )

                def func(data2d, labels2d):
                    datares, opts = c_filter.charge_artifact_FD_filter_downup_av_prevlines3_2d(data2d, labels2d)
                    return datares

                cCalculate = cSliceBySliceFunctHelper(self.viewer,func, curData, curLabels)
            
            if not cCalculate is None:

                #TODO: USe the napari Activity progress bar

                estimate_iters = cCalculate.get_len()
                print(f"estimate_iters:{estimate_iters}")
                self.iter=0
                def calbackfn():
                    print(f"slice iter:{self.iter} / {estimate_iters}")
                    self.iter+=1

                data_filt_whole = cCalculate.run(calbackfn)

                print("Calculation complete")

                if not data_filt_whole is None:
                    self.viewer.add_image(data_filt_whole, name="FFT-dir filter")


        # if self.chkPreview.isChecked():
        #     #Check if there is valid data for creating the preview
            
        #     if not self.curDataSlice is None:
        #         #Debug
        #         #print(f"curData.shape = {self.curData.shape}")
        #         self.labelSliceDims.setText(str(self.curDataSlice.shape))

        #         freq_low = self.freqselector_low.value()
        #         freq_high = self.freqselector_high.value()

        #         #print(f"freq_low : {freq_low} , freq_high: {freq_high}")

        #         #Calculate the filtered data
        #         dataPreview = filters.fft_bandpass_filter_dirx_2D(self.curDataSlice,freq_low, freq_high)

        #         if self.myImageLayer is None or not self.myImageLayer in self.viewer.layers:
        #             #self.myImageLayer = self.viewer.add_image(dataPreview, name="preview")
        #             #self.myImageLayer._keep_auto_contrast=True #This helps show image with auto contrast however associated button does not show it is pressed
        #             #self.myImageLayer.reset_contrast_limits() #Try, doesnt work
        #             self.myImageLayer = self.viewer.add_image(dataPreview, name="preview", contrast_limits= (self.curDataSlice.min(), self.curDataSlice.max()))

        #         else:
        #             #Dont create a new layer, just modify it
        #             self.myImageLayer.data = dataPreview 
        #             #self.myImageLayer.refresh()

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

        data0 = self.w_setselect.get_active_image_selected_data_slice()

        curDataSlice=data0
        if not curDataSlice is None:
            #Calculate
            overlay = measure_charging.generate_heatmap(curDataSlice, tilesize)
            self.viewer.add_image(overlay, name="Charging artefacts", opacity=0.3, colormap="viridis")


    def btnSACalculate_onclick(self):
        print("btnSACalculate_onclick()")

        res=None

        data3d = self.w_setselect.get_active_image_selected_data()

        if data3d.ndim==3:

            #Estimate number of iterations for progress bar
            nslices = data3d.shape[0]
            niterations = 2*(nslices+1)

            pbr=progress(total=niterations,desc="Alignment progress")

            def callbkfn():
                pbr.update(1)
                pbr.refresh()

            #Show or not show Activity dialog?
            #activ_dialog=self.viewer.window._qt_viewer.window()._activity_dialog #Warning this will be unavailable in the future
            #activ_dialog.show()

            res= slice_alignment.align_stack(data3d, slice_alignment.ALIGNMENT_METHOD_DEFAULT,callbkfn)
            
            pbr.close()
            #activ_dialog.hide()

            if not res is None:
                self.viewer.add_image(res, name="stack aligned")

    def btnQuollCalcFRC_onclick(self):
        print("btnQuollCalcFRC_onclick()")

        with progress(total=2,desc="FRC calculation in progress") as pbr:

            data2d = self.w_setselect.get_active_image_selected_data_slice()

            pbr.update(1)

            res = quoll.getFRC(data2d)

            pbr.update(1)

            print (f"FRC resolution : {res}")

            #print result on a widget
            self.lblFRCWholeImg.setText( f"{res:.4e}")
        return


    def btnQuollCalcFRCTiled_onclick(self):
        print("btnQuollCalcFRCTiled_onclick()")

        #Gathers tile size parameter
        tile_size = self.spbQuolTileSize.value()

        #Runs the FRC calculation for tiles
        with progress(total=2,desc="FRC tiled calculation in progress") as pbr:

            data2d = self.w_setselect.get_active_image_selected_data_slice()

            pbr.update(1)

            res = quoll.getTiledFRC(data2d, tilesize=tile_size)

            pbr.update(1)

            print (f"FRC tiled result : {res}")

            if not res is None:
                dfres, heatmap = res

                print(f"FRC result summary:{dfres.describe().to_dict()}")
                #self.frctable = widget_PandasDFTable(self.viewer, dfres)
                #self.quollTableContainerWidget.layout().addWidget(self.frctable)

                #self.frctable = widget_PandasDFTable(self.quollTableContainerWidget, dfres)
                #TODO remove all widgets if any
                #self.quollTableContainerWidgetLayout.addWidget(widget_PandasDFTable(self.viewer, dfres))
                #self.frctable = widget_PandasDFTable(self.viewer, dfres)
                self.quollTable.setDataframe(dfres) #hopefully it will update by itself
                
                #TODO: display a heatmap
                self.viewer.add_image(heatmap,name="FRC heatmap",
                    colormap="viridis",
                    opacity=0.3
                )
        return

                


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






class widget_ImageSetCurrentSelect(QtWidgets.QWidget):
    '''
    A QtWidget that helps select Image data for further processing

    Needs to be initialised with a napari Viewer object

    '''

    def __init__(self, napari_viewer: 'napari.viewer.Viewer', header=""):
        super().__init__()

        self.viewer = napari_viewer #Reference to the viewer that will be needed later
        #self.viewer = napari.current_viewer() #Doesnt work, returns None

        layout0 = QtWidgets.QVBoxLayout()
        self.setLayout(layout0)

        group_box0 = QtWidgets.QGroupBox(header)
        #self.t0_vbox0.addWidget(self.group_box0)
        layout0.addWidget(group_box0)

        gridlayout = QtWidgets.QGridLayout()
        group_box0.setLayout(gridlayout)

        self.b_set = QtWidgets.QPushButton("Set currently selected")
        gridlayout.addWidget(self.b_set,0,0)
        self.b_set.clicked.connect(self.b_set_on_click) #Signal
  
        gridlayout.addWidget(QtWidgets.QLabel("name"),1,0)
        self.l_name_v = QtWidgets.QLabel("")
        gridlayout.addWidget(self.l_name_v,1,1)

        gridlayout.addWidget(QtWidgets.QLabel("shape"),2,0)
        self.l_shape_v = QtWidgets.QLabel("")
        gridlayout.addWidget(self.l_shape_v,2,1)

        #layout.addWidget(QtWidgets.QLabel("slicing axis"),3,0)
        self.l_slicingaxis_v = QtWidgets.QLabel("")
        # layout.addWidget(self.l_slicingaxis_v,3,1)
        # #self.l_slicingaxis_v.hide() #Hide
        # self.l_slicingaxis_v.setHidden(True) #Hide

        #layout.addWidget(QtWidgets.QLabel("slice shape"),4,0)
        self.l_sliceshape_v = QtWidgets.QLabel("")
        # layout.addWidget(self.l_sliceshape_v,4,1)
        # self.l_sliceshape_v.hide() #Hide

        #layout.addWidget(QtWidgets.QLabel("nslice"),5,0)
        self.l_nslice_v = QtWidgets.QLabel("")
        # layout.addWidget(self.l_nslice_v,5,1)
        # self.l_nslice_v.hide() #Hide

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
        print(f"type(nap_image) : {type(nap_image) }")
        if isinstance(nap_image, napari.layers.Image) or isinstance(nap_image, napari.layers.Labels):
            self.l_name_v.setText(nap_image.name)
            print("Current active layer name is: ", nap_image.name)
            self.curImage = nap_image
            self.l_shape_v.setText(str(self.curImage.data.shape))

    def get_active_image(self):
        active0=None
        if not self.viewer is None:
            active0 = self.viewer.layers.selection.active #check if any layer is active
            self.set_default_napari_image(active0)
        else:
            print("get_active_image(): no napari.viewer, cannot get active image")

        return active0

    def get_active_image_selected_data(self) ->ImageData:
        print("get_active_image_selected_data()")
        if self.curImage is None:
            self.b_set_on_click() #Does the same thing as clicking the button to set default image data
            #return self.curImage.data
        
        return self.curImage.data

    def get_active_image_selected_data_slice(self):
        print("get_active_image_selected_data_slice()")

        curDataSlice=None

        curData = self.get_active_image_selected_data()
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
 


    # #TODO: Create function that returns iterator
    # class iteratorSliceBySlice:
    #     '''
    #     Iterator that gets data slice by slice
        
    #     It is implemented using the alternative iterator protocol
    #     that uses __getitem__ and __len__'''
    #     def __init__(self,data, viewer):
    #         if data.ndim==3:
    #             #gets orientation of viewing. If viewer is in 3D it sets to z-axis order=0
    #             self.order=0
    #             if self.viewer.dims.ndisplay == 2:
    #                 dims0 = self.viewer.dims
    #                 self.order = dims0.order[0] #Gets the axis (or plane), that the Viewer is currently viewing
                

    #             self.i=0

    #             self.data_oriented = None

    #             if self.order == 0: #z
    #                 self.data_oriented=data
    #                 self.l_slicingaxis_v.setText('z')
    #             elif self.order == 1: #y
    #                 self.data_oriented = np.transpose(data,(1,2,0))
    #                 self.l_slicingaxis_v.setText("y")
    #             else: #x
    #                 self.data_oriented = np.transpose(data,(2,0,1))
    #                 self.l_slicingaxis_v.setText("x")

                
    #         else:
    #             return None

    #     def __iter__(self, data, viewer):
            

    #         return self #Needed

    #     def __next__(self):

class cSliceBySliceFunctHelper():

    def __init__(self,
            napari_viewer: 'napari.viewer.Viewer',
            func,
            data, 
            data_opt=None 
        ):
        '''
            data
            func to apply to slices. Must be like func(data2D) -> datares2D

            Some functions require two datasets, hence the optional data3d_opt
            In this case func(data2D, data3d_opt) -> datares2D
        '''

        #viewer = napari.current_viewer()
        viewer = napari_viewer

        self.datadim = data.ndim
        self.data_oriented = None
        self.data_opt_oriented = None

        if data.ndim==3:
            #gets orientation of viewing. If viewer is in 3D it sets to z-axis order=0
            self.order=0
            if viewer.dims.ndisplay == 2:
                dims0 = viewer.dims
                self.order = dims0.order[0] #Gets the axis (or plane), that the Viewer is currently viewing

            self.has_opt = False
            if not data_opt is None:
                self.has_opt = True

            if self.order == 0: #z
                self.data_oriented=data
                #self.l_slicingaxis_v.setText('z')
                if self.has_opt:
                    self.data_opt_oriented=data_opt

            elif self.order == 1: #y
                self.data_oriented = np.transpose(data,(1,2,0))
                #self.l_slicingaxis_v.setText("y")
                if self.has_opt:
                    self.data_opt_oriented= np.transpose(data_opt,(1,2,0))
            else: #x
                self.data_oriented = np.transpose(data,(2,0,1))
                #self.l_slicingaxis_v.setText("x")
                if self.has_opt:
                    self.data_opt_oriented= np.transpose(data_opt,(2,0,1))

            self.n= self.data_oriented.shape[0]

            self.func = func

        else:
            print("data is not 3D. Function will be applied to the single slice")
            self.data_oriented = data
            self.data_opt_oriented = data_opt
            self.n=1
    
    def run(self, callbackfn=None):
        if self.datadim==3:
            datares_shaperestored=None
            
            if not self.data_oriented is None:
                
                #Check result type with first slice
                if not self.has_opt:
                    slice0res = self.func(self.data_oriented[0,:,:])
                else:
                    slice0res = self.func(self.data_oriented[0,:,:], self.data_opt_oriented[0,:,:])

                #Set result to have the same type as first layer result
                datares = np.zeros_like(self.data_oriented, dtype=slice0res.dtype)

                datares[0,:,:]= slice0res

                #callback
                if not callbackfn is None:
                    callbackfn()

                #Do for the rest of slices
                for i in range(1, self.n):
                    #dataslice = self.data_oriented[i,:,:]

                    #Apply function
                    if not self.has_opt:
                        datasliceres = self.func(self.data_oriented[i,:,:])
                    else:
                        datasliceres = self.func(self.data_oriented[i,:,:],self.data_opt_oriented[i,:,:])
                    
                    datares[i,:,:] = datasliceres

                    #callback
                    if not callbackfn is None:
                        callbackfn()
                    
                #Upon completion restore the shape of the data
                if self.order == 0: #z
                    datares_shaperestored = datares
                elif self.order == 1: #y
                    datares_shaperestored = np.transpose(datares,(2,0,1))
                else: #x
                    datares_shaperestored = np.transpose(datares,(1,2,0))

            return datares_shaperestored
        
        elif self.datadim==2:
            data2dres = None
             #Apply function
            if not self.has_opt:
                data2dres = self.func(self.data_oriented)
            else:
                data2dres = self.func(self.data_oriented,self.data_opt_oriented)
            return data2dres
        
        return None
    
    def get_len(self):
        #Returns number of iterations expected for the calculation

        niter=0
        if not self.data_oriented is None:
            niter= self.n

        return niter


#NOT WORKING
# Try to write a widget for chaffer using magicgui

# # In Redlionfish package the decorator @magicfactory is used instead
# #@magicgui(
# @magic_factory(
#     imgd1={'label': 'Data'}, 
#     imgd2={'label': 'Labels'},
#     call_button="execute",
#     layout="vertical",
#     bThisSlice={'label':'This slice'}
#     )
# def chaffer_magicgui_widget(
#         imgd1: ImageData,
#         imgd2: ImageData,
#         bThisSlice=True,
#         nlinesaverage=20,
#         data_fit_max_length_edge_px = 700,
#         data_fit_min_length_px = 50
#     ) -> LayerDataTuple: #Result is a LayerDataTuple, like (data, {dict_properties})

#     datares = None

#     c_filter = chafer.cls_charge_artifact_suppression_filter(nlinesaverage,data_fit_max_length_edge_px, data_fit_min_length_px)

#     if bThisSlice:
#         res0, opts = c_filter.charge_artifact_FD_filter_downup_av_prevlines3_2d(imgd1, imgd2)
#         ret = ( res0 , { 'name':'chafer result of slice'})
#     else:
#         res0, opts = c_filter.charge_artifact_FD_filter_downup_av_prevlines3_2d(imgd1, imgd2)
#         ret = ( res0 , { 'name':'chafer result'})

#     return ret


class widget_PandasDFTable(QtWidgets.QTableWidget):
    '''
    A simple table that displays a pandas dataframe without interaction

    '''
    def __init__(self, napari_viewer: 'napari.viewer.Viewer'):
        super().__init__()

    def setDataframe(self,dataframe):
        self.clearContents()

        self.dataframe = dataframe

        #Fill table with content from dataframe
        rowcount = dataframe.shape[0]
        columncount = dataframe.shape[1]

        self.setRowCount(rowcount+1)
        self.setColumnCount(columncount)

        #Header
        for i,colnames in enumerate(dataframe.keys()):
            self.setItem(0,i, QtWidgets.QTableWidgetItem(str(colnames)))
        
        #data
        for i, colnames in enumerate(dataframe.keys()):
            for j, value in enumerate(dataframe.get(colnames)):
                self.setItem(j, i, QtWidgets.QTableWidgetItem(str(value)))
