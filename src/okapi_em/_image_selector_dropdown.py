from qtpy import QtWidgets
from qtpy.QtCore import Qt

import napari

class image_selector_dropdown_widget(QtWidgets.QComboBox):

    # Dropdown combobox that contains images

    def __init__(self, viewer, *args, **kwargs):
        super(QtWidgets.QComboBox, self).__init__(*args, **kwargs)


