#Basic publications set up.

#Run manually with something like
#sudo puppet apply --verbose ./papers-setup.pp

#Olly Butters
#19/6/17

#todo
##Do group membership better - users are hardcoded in here.

node 'pubs-alspac-p0' {

	#Root path everything done wrt
	$root_path = '/var/local/projects/puma/'

	#############################################
	#Folder tree first
	#############################################

	#Top level folder where everything goes.
	file {'root_dir':
		path   => $root_path,
		ensure => directory,
	}

	#Apache
	file {'apache_dir':
		path    => "${root_path}apache",
    ensure  => directory,
		require => File['root_dir']
  }

	#Apache conf
	file {'apach_conf_dir':
		path    => "${root_path}apache/conf",
    ensure  => directory,
		owner   => 'root',
		group   => 'd2k_devs',
		mode    => '660',
		seluser => 'system_u',
		selrole => 'object_r',
		seltype => 'httpd_config_t',
		recurse => true,
		require => [Group['d2k_devs'], File['apache_dir']],
	}

	#Apache logs
	file {'apache_logs_dir':
		path    => "${root_path}apache/logs",
		ensure  => directory,
		owner   => 'root',
		group   => 'd2k_devs',
		mode    => '640',
		seluser => 'system_u',
		selrole => 'object_r',
		seltype => 'httpd_log_t',
		recurse => true,
		require => [Group['d2k_devs'], File['apache_dir']],
	}

	#bin
	file {'bin_dir':
		path    => "${root_path}bin",
    ensure  => directory,
		owner   => 'root',
		group   => 'd2k_devs',
		mode    => '775',
		require => [Group['d2k_devs'], File['root_dir']],
  }

	#www
	file {'www_dir':
		path    => "${root_path}www",
    ensure  => directory,
		owner   => 'root',
		group   => 'd2k_devs',
		seluser => 'system_u',
		selrole => 'object_r',
		seltype => 'httpd_sys_content_t',
		recurse => true,
		require => [Group['d2k_devs'], File['root_dir']],
  }

	#www - puma
	#The logic here is that root owns the files, so effectively the owner never looks at them.
	#The d2k_devs group has everyone in it that would ever edit any www files, and everyone else
	#i.e. apache has r-x access on all the files.
	file {'www_puma_dir':
		path    => "${root_path}www/puma",
    ensure  => directory,
		owner   => 'root',
		group   => 'd2k_devs',
		mode    => '775',
		seluser => 'system_u',
		selrole => 'object_r',
		seltype => 'httpd_sys_content_t',
		recurse => true,
		require => [File['www_dir'],Group['d2k_devs']],
  }


	##############################################
	#Packages
	##############################################

	#Apache
	package {'httpd':
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


	################################################
	#Users and groups
	################################################

	#Make a group for the developers
	group {'d2k_devs':
		ensure  => 'present',
	}

	user {'nob22':
		groups  => 'd2k_devs',
		require => Group['d2k_devs'],
	}
}
