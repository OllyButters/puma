#Basic publications set up.

#Run manually with something like
#sudo puppet apply --verbose ./papers-setup.pp

#Olly Butters
#10/12/16

#todo
##Add a variable to manage the path so it becomes a little more portable.
##Do group membership better - users are hardcoded in here.

node 'pubs-alspac-p0' {

	#############################################
	#Folder tree first
	#############################################

	#Top level folder where everything goes.
	file {'root_dir':
		path   => '/var/local/projects/alspac',
		ensure => directory,
	}

	#Apache
	file {'apache_dir':
		path    => '/var/local/projects/alspac/apache',
                ensure  => directory,
		require => File['root_dir']
        }

	#Apache conf
	file {'apach_conf_dir':
		path    => '/var/local/projects/alspac/apache/conf',
                ensure  => directory,
		seluser => 'system_u',
		selrole => 'object_r',
		seltype => 'httpd_config_t',
		recurse => true,
		require => File['apache_dir'],
        }

	#Apache logs
	file {'apache_logs_dir':
		path    => '/var/local/projects/alspac/apache/logs',
		ensure  => directory,
		owner   => 'root',
		group   => 'alspac_devs',
		mode    => '640',
		seluser => 'system_u',
		selrole => 'object_r',
		seltype => 'httpd_log_t',
		recurse => true,
		require => File['apache_dir'],
	}

	#bin
	file {'bin_dir':
		path    => '/var/local/projects/alspac/bin',
                ensure  => directory,
		owner   => 'root',
		group   => 'alspac_devs',
		mode    => '775',
		require => [Group['alspac_devs'], File['root_dir']],
        }

	#keys
	file {'keys_dir':
		path    => '/var/local/projects/alspac/keys',
                ensure  => directory,
		require => File['root_dir']
        }

	#nagios
	file {'nagios_dir':
		path    => '/var/local/projects/alspac/nagios',
                ensure  => directory,
		require => File['root_dir'],
        }

	#var
	file {'var_dir':
		path    => '/var/local/projects/alspac/var',
                ensure  => directory,
		require => File['root_dir']
        }

	#www
	file {'www_dir':
		path    => '/var/local/projects/alspac/www',
                ensure  => directory,
		owner   => 'root',
		group   => 'alspac_devs',
		seluser => 'system_u',
		selrole => 'object_r',
		seltype => 'httpd_sys_content_t',
		recurse => true,
		require => [Group['alspac_devs'], File['root_dir']],
        }

	#www - pubs
	#The logic here is that root owns the files, so effectively the owner never looks at them.
	#The alspac_devs group has everyone in it that would ever edit any www files, and everyone else
	#i.e. apache has r-x access on all the files. 
	file {'www_pubs_dir':
		path    => '/var/local/projects/alspac/www/pubs',
                ensure  => directory,
		owner   => 'root',
		group   => 'alspac_devs',
		mode    => '775',
		seluser => 'system_u',
		selrole => 'object_r',
		seltype => 'httpd_sys_content_t',
		recurse => true,
		require => [File['www_dir'],Group['alspac_devs']],
        }


	##############################################
	#Packages
	##############################################

	#Apache
	package {'httpd':
		ensure => 'installed',
	}

	#mod ssl
	package {'mod_ssl':
		ensure => 'installed',
	}

	#Nagios plugins
	package {'nagios-plugins-all':
		ensure => 'installed',
	}

	#Make sure apache is running	
	service {'httpd':
		ensure => 'running',
		enable => 'true',
	}

	#GCC needed to build some python libraries
	package {'python-devel':
    		ensure => 'installed',
	}

	#Python libraries
	package {'python2-pip':
    		ensure  => 'installed',
		require => Package['python-devel'],
	}

	#NOTE: there is a bug in pip in this version of CENTOS where it looks
	#for pip in /usr/bin/pip-python, but that doesnt exist. Put in a soft
	#link to where pip actually is /usr/bin/pip. Ive not done this as a 
	#puppet thing as I hope it will get fixed.

	#Python libraries

	#biopython
	package {'biopython':
    		ensure => 'installed',
    		provider => 'pip',
		require => Package['python2-pip'],		
	}

	#htmlentities
	package {'htmlentities':
    		ensure => 'installed',
    		provider => 'pip',
		require => Package['python2-pip'],		
	}


	#jsonpath-rw
	package {'jsonpath-rw':
    		ensure => 'installed',
    		provider => 'pip',
		require => Package['python2-pip'],
	}

	#requests
	package {'requests':
    		ensure => 'installed',
    		provider => 'pip',
		require => Package['python2-pip'],		
	}

	#pyzotero
	package {'pyzotero':
    		ensure => 'installed',
    		provider => 'pip',
		require => Package['python2-pip'],		
	}

	#unicodecsv
	package {'unicodecsv':
    		ensure => 'installed',
    		provider => 'pip',
		require => Package['python2-pip'],		
	}


	#################################################
	#Misc SELinux settings
	#################################################
	selboolean {'httpd_can_network_connect':
		persistent => true,
		value      => on,
	}


	#################################################
	#CRON jobs
	#################################################
	#nagios check disk usage
	cron { nagios_check_disk:
		command => '/var/local/projects/alspac/nagios/check_disk.sh > /dev/null 2>&1',
		user    => root,
		minute  => '*/5',
	}

	#nagios check load
	cron { nagios_check_load:
		command => '/var/local/projects/alspac/nagios/check_load.sh > /dev/null 2>&1',
		user    => root,
		minute  => '*/5',
	}
	
	#nagios check number of OS users logged in
	cron { nagios_check_users:
		command => '/var/local/projects/alspac/nagios/check_users.sh > /dev/null 2>&1',
		user    => root,
		minute  => '*/5',
	}

	#nagios check total number of processes running
	cron { nagios_check_total_procs:
		command => '/var/local/projects/alspac/nagios/check_total_procs.sh > /dev/null 2>&1',
		user    => root,
		minute  => '*/5',
	}

	#nagios check ntp running
	cron { nagios_check_ntp:
		command => '/var/local/projects/alspac/nagios/check_ntp.sh > /dev/null 2>&1',
		user    => root,
		minute  => '*/5',
	}


	################################################
	#Users and groups
	################################################

	#Make a group for the developers
	group {'alspac_devs':
		ensure  => 'present',
	}

	user {'epxol':
		groups  => 'alspac_devs',
		require => Group['alspac_devs'],
	}

	user {'ob13747':
		groups  => 'alspac_devs',
		require => Group['alspac_devs'],
	}

	user {'epnsw':
		groups  => 'alspac_devs',
		require => Group['alspac_devs'],
	}

	user {'edzhg':
		groups  => 'alspac_devs',
		require => Group['alspac_devs'],
	}
}
