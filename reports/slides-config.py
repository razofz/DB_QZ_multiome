# `my_slide_config.py` 
c.TagRemovePreprocessor.remove_input_tags.append("to_remove")
c.SlidesExporter.reveal_theme="serif"
c.SlidesExporter.reveal_number = 'c/t'
    # slide number format (e.g. ‘c/t’). Choose from: ‘c’: current, ‘t’: total, ‘h’: horizontal, ‘v’: vertical
c.SlidesExporter.reveal_scroll = True
    # If True, enable scrolling within each slide
c.SlidesExporter.reveal_transition = 'slide'
    # The list of transitions that ships by default with reveal.js are: none, fade, slide, convex, concave and zoom.
c.SlidesExporter.theme = 'light'
print(c)
print(dir(c))
print(type(c))

# the following does the equivalent of --to slides and --post serve, see here: https://github.com/jupyter/nbconvert/blob/master/setup.py#L220 + https://github.com/jupyter/nbconvert/blob/master/nbconvert/nbconvertapp.py#L482 and here https://github.com/jupyter/nbconvert/blob/master/nbconvert/nbconvertapp.py#L238
app_settings = {
    # "postprocessor_class": "nbconvert.postprocessors.serve.ServePostProcessor",
    "export_format": "slides"
}
c.NbConvertApp.update(app_settings)

# the following does the equivalent of --no-prompt, see here: https://github.com/jupyter/nbconvert/blob/master/nbconvert/nbconvertapp.py#L109-L114
exporter_settings = {
    'exclude_input_prompt' : True,
    'exclude_input' : True,
    'exclude_output_prompt' : True,
}
c.TemplateExporter.update(exporter_settings)
