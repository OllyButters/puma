import os
import shutil
from config import config

# Copy specified reports to the html directory
def copy_reports():
    print('Copying reports to html directory.')

    if config.WEB_PAGE_REPORTS:
        print('Copying ', config.WEB_PAGE_REPORTS, ' to html directory.')

        if os.path.exists(config.data_dir + '/' + config.WEB_PAGE_REPORTS):

            shutil.copy(config.data_dir + '/' + config.WEB_PAGE_REPORTS, config.html_dir + '/reports/')
