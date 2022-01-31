import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

dfctpc = pd.DataFrame(celltags_per_cell, columns=['cell_count'])
dfctpc['celltag_count'] = dfctpc.index

fig, ax = plt.subplots()
ax.bar(dfctpc.celltag_count, height=dfctpc.cell_count)
ax.set_title(f'Distribution of cells per CellTag for {sample}')

ax.set_ylabel('Cells with this number of CellTags')
ax.set_xlabel('CellTags per cell')

plt.savefig('cellsperCellTag.png')



dfcpct = pd.DataFrame(cells_per_celltag, columns=['cell_count'])
dfcpct['celltag_count'] = dfcpct.index

fig, ax = plt.subplots()
ax.bar(dfcpct.celltag_count, height=dfcpct.cell_count)
ax.set_title(f'Distribution of CellTags per cell for {sample}')

ax.set_ylabel('CellTags with this number of cells')
ax.set_xlabel('Cells per CellTags ')

plt.savefig('CellTagspercell.png')


# def autolabel(rects):
#     """
#     Attach a text label above each bar displaying its height
#     """
#     for rect in rects:
#         height = rect.get_height()
#         ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
#                 '%d' % int(height),
#                 ha='center', va='bottom')


# autolabel(rects)
