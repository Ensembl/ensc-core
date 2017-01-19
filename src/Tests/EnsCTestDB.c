/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


=pod

=head1 NAME - OtterTestDB

=head1 SYNOPSIS

    # Add test dir to lib search path
    use lib 't';
    
    use OtterTestDB;
    
    my $otter_test = OtterTestDB->new();
    # OR: 
    #   my $otter_test = OtterTestDB->new({host=>'foo'});
    # OR:
    #   my $otter_test = OtterTestDB->new('myconf.dat');
    
    # Load some data into the db
    $otter_test->do_sql_file("some_data.sql");
    
    # Get an Otter db object for the test db
    my $otterdb = $otter_test->get_DBSQL_Obj;

=head1 DESCRIPTION

This is a module used just by the Otter test
suite to create a test database for a particular
test.  Creating a new object creates a database
with a name such that it should never clash with
other users testing on the same server.  The
database is destroyed when the object goes out of
scope.

The settings, such as the server host and port,
are found in the file B<OtterTestDB.conf>.  See
B<OtterTestDB.conf.example> for an example.


Firstly, the file OtterTestDB.conf will be read. If this fails, some
(hopefully reasonable) defaults are taken. Secondly, if an argument is
given, this will be used to get arguments from. They will supplement or
override those of OtterTestDB.conf. If the argument is a filename, that file
will be read, and its contents will override/merge with those of
OtterTestDB.conf. If the argument is a hash, its contents will likewise
override/merge with those read from OtterTestDB.conf

=head1 METHODS

=cut

package OtterTestDB;

use Bio::Otter::DBSQL::DBAdaptor;
use vars qw(@ISA);
use strict;
use Sys::Hostname 'hostname';
use DBI;
use Carp;

#Package variable for unique database name
my $counter=0;

my $static_test_dbh = undef;
my $static_test_dbname = "";

{
    # This is a list of possible entries in the config
    # file "OtterTestDB.conf" or in the hash being used.
    my %known_field = map {$_, 1} qw(
        driver
        host
        user
        port
        pass
        schema_sql
        module
        );

    ### Firstly, the file OtterTestDB.conf will be read. If this fails, some
    ### hopefully reasonable defaults are taken. Secondly, if an argument
    ### is given, this will be used to get arguments from. They will
    ### supplement or override those of OtterTestDB.conf. If the argument
    ### is a filename, that file will be read, and its contents will
    ### override/merge with those of OtterTestDB.conf. If the argument is a
    ### hash, its contents will likewise override/merge with those read
    ### from OtterTestDB.conf
    sub new {
        my( $pkg, $arg ) = @_;

        $counter++;
        my $conf_file='OtterTestDB.conf';
        my $fallback_defaults = {
		'driver'        => 'mysql',
		'host'          => 'localhost',
		'user'          => 'root',
		'port'          => '3306',
		'pass'      => undef,
		'schema_sql'    => ['../sql/table.sql'],
		'module'        => 'Bio::Otter::DBSQL::DBAdaptor'
		};

        my $self =undef;
	$self = do $conf_file || $fallback_defaults;
	if ($arg) {
	    if  (ref $arg eq 'HASH' ) {  # a hash ref
                foreach my $key (keys %$arg) {
		    $self->{$key} = $arg->{$key};
		}
	    }
	    elsif (-f $arg )  { # a file name
		$self = do $arg;
	    } else {
		confess "expected a hash ref or existing file";
	    }
	}
        
        foreach my $f (keys %$self) {
            confess "Unknown config field: '$f'" unless $known_field{$f};
        }
        bless $self, $pkg;
        $self->create_db;
	
        return $self;
    }
}

sub driver {
    my( $self, $value ) = @_;
    
    if ($value) {
        $self->{'driver'} = $value;
    }
    return $self->{'driver'} || confess "driver not set";
}

sub host {
    my( $self, $value ) = @_;
    
    if ($value) {
        $self->{'host'} = $value;
    }
    return $self->{'host'} || confess "host not set";
}

sub user {
    my( $self, $value ) = @_;
    
    if ($value) {
        $self->{'user'} = $value;
    }
    return $self->{'user'} || confess "user not set";
}

sub port {
    my( $self, $value ) = @_;
    
    if ($value) {
        $self->{'port'} = $value;
    }
    return $self->{'port'};
}

sub pass {
    my( $self, $value ) = @_;
    
    if ($value) {
        $self->{'pass'} = $value;
    }
    return $self->{'pass'};
}

sub schema_sql {
    my( $self, $value ) = @_;
    
    if ($value) {
        push(@{$self->{'schema_sql'}}, $value);
    }
    return $self->{'schema_sql'} || confess "schema_sql not set";
}

sub dbname {
    my( $self ) = @_;

    $self->{'_dbname'} ||= $self->_create_db_name();
    return $self->{'_dbname'};
}

# convenience method: by calling it, you get the name of the database,
# which  you can cut-n-paste into another window for doing some mysql
# stuff interactively
sub pause {
    my ($self) = @_;
    my $db = $self->{'_dbname'};
    print STDERR "pausing to inspect database; name of database is:  $db\n";
    print STDERR "press ^D to continue\n";
    `cat `;
}

sub module {
    my ($self, $value) = @_;
    $self->{'module'} = $value if ($value);
    return $self->{'module'};
}

sub _create_db_name {
    my( $self ) = @_;

    my $host = hostname();
    my $db_name = "_test_db_${host}_$$".$counter;
    $db_name =~ s{\W}{_}g;
    return $db_name;
}

sub create_db {
    my( $self ) = @_;
    
    ### FIXME: not portable between different drivers
    my $locator = 'DBI:'. $self->driver .':host='.$self->host;
    if( defined $self->port() ) {
      $locator .= ";port=".$self->port();
    }
    my $db = DBI->connect(
        $locator, $self->user, $self->pass, {RaiseError => 1}
        ) or confess "Can't connect to server";
    my $db_name = $self->dbname;
    $db->do("CREATE DATABASE $db_name");
    $db->disconnect;
    
    
    $self->do_sql_file(@{$self->schema_sql});
}

sub db_handle {
    my( $self ) = @_;
    
    unless ($self->{'_db_handle'}) {
        $self->{'_db_handle'} = DBI->connect(
            $self->test_locator, $self->user, $self->pass, {RaiseError => 1}
            ) or confess "Can't connect to server";
    }
    return $self->{'_db_handle'};
}

sub test_locator {
    my( $self ) = @_;
    
    my $locator = 'dbi:'. $self->driver .':database='. $self->dbname;
    foreach my $meth (qw{ host port }) {
        if (my $value = $self->$meth()) {
            $locator .= ";$meth=$value";
        }
    }
    return $locator;
}

sub otter_locator {
    my( $self) = @_;
    
    my $module = ($self->module() || 'Bio::Otter::DBSQL::DBAdaptor');
    my $locator = '';
    foreach my $meth (qw{ host port dbname user pass }) {
        my $value = $self->$meth();
	next unless defined $value;
        $locator .= ';' if $locator;
        $locator .= "$meth=$value";
    }
    $locator .= ";perlonlyfeatures=1";
   
    return "$module/$locator";
}

# return the database handle:
sub get_DBSQL_Obj {
    my( $self ) = @_;
    
    my $locator = $self->otter_locator();

    my $db =  new Bio::Otter::DBSQL::DBAdaptor(-host => $self->host,
	-user => $self->user,
	-pass => $self->pass,
	-port => $self->port,
	-dbname => $self->dbname);
}

sub do_sql_file {
    my( $self, @files ) = @_;
    local *SQL;
    my $i = 0;
    my $dbh = $self->db_handle;

    my $comment_strip_warned=0;

    foreach my $file (@files)
    {
        my $sql = '';
        open SQL, $file or die "Can't read SQL file '$file' : $!";
        while (<SQL>) {
            # careful with stripping out comments; quoted text
            # (e.g. aligments) may contain them. Just warn (once) and ignore
            if (    /'[^']*#[^']*'/ 
                 || /'[^']*--[^']*'/ ) {
                     if ( $comment_strip_warned++ ) { 
                         # already warned
                     } else {
                         warn "#################################\n".
                           warn "# found comment strings inside quoted string; not stripping, too complicated: $_\n";
                         warn "# (continuing, assuming all these they are simply valid quoted strings)\n";
                         warn "#################################\n";
                     }
                 } else {
                s/(#|--).*//;       # Remove comments
            }
            next unless /\S/;   # Skip lines which are all space
            $sql .= $_;
            $sql .= ' ';
        }
        close SQL;
        
	#Modified split statement, only semicolumns before end of line,
	#so we can have them inside a string in the statement
	#\s*\n, takes in account the case when there is space before the new line
        foreach my $s (grep /\S/, split /;[ \t]*\n/, $sql) {
	    $s =~ s/\;\s*$//g;
            $self->validate_sql($s);
            $dbh->do($s);
            $i++
        }
    }
    return $i;
}                                       # do_sql_file

sub validate_sql {
    my ($self, $statement) = @_;
    if ($statement =~ /insert/i)
    {
        $statement =~ s/\n/ /g; #remove newlines
        die ("INSERT should use explicit column names (-c switch in mysqldump)\n$statement\n")
            unless ($statement =~ /insert.+into.*\(.+\).+values.*\(.+\)/i);
    }
}

sub dnadb {
  my ($self,$dnadb) = @_;

  if (defined($dnadb)) {
     $self->{_dnadb} = $dnadb;
  }
  return $self->{_dnadb};
}


sub DESTROY {
    my( $self ) = @_;

    if (my $dbh = $self->db_handle) {
        my $db_name = $self->dbname;
        $dbh->do("DROP DATABASE $db_name");
        $dbh->disconnect;
    }
}


=head2 test_db

  Arg 1     : txt $database_name
              is not autogenerated, must be removed explicitly
  Function  : returns the test-db filled with approx 2MB of data in the
              current schema. Done once on first call. Uploads test-db
              subdir into normal test database.
  Returntype: Bio::Otter::DBSQL::DBAdaptor
  Exceptions: none
  Caller    : test modules

=cut

sub test_db {
   my $self = shift;
   my $dbname = shift;

   if( ! defined $static_test_dbh ) {
     # dbhandle generation ...
     my $locator = 'DBI:'. $self->driver .':host='.$self->host;
     if( defined $self->port() ) {
       $locator .= ";port=".$self->port();
     }

     my $db = DBI->connect(
			   $locator, $self->user, $self->pass, {RaiseError => 1}
			  ) or confess "Can't connect to server";
     my $sth = $db->prepare( "show databases" );
     my $dbexists = 0;
     $sth->execute();

     while( my $aref = $sth->fetchrow_arrayref ) {
       if( $aref->[0] eq $dbname ) {
	 # db exists
	 $dbexists = 1;
	 last;
       }
     }

     # check if db exists, then dont upload
     if( $dbexists ) {

       $db->do( "use $dbname" );
       $static_test_dbh = $db;
       $static_test_dbname = $dbname;

     } else {
       # create and upload database
       $db->do( "create database $dbname"  );
       $db->do( "use $dbname" );

       my $dir = __FILE__;

       $dir =~ s/[^\/]*$//g;
       my @files = glob( $dir."test-db/*.sql" );
       my $tmp_db = $self->db_handle();

       # create the schema
       $self->{_db_handle} = $db ;
       $self->do_sql_file( @files );
       $self->{_db_handle}= $tmp_db;

       # import the data
       for my $file ( @files ) {
	 my $datafile = $file;
	 $datafile =~ s/\.sql$/.txt/g;
	 my ( $tablename ) = $file =~ /([^\/]*)\.sql/;
	 $db->do( "load data local infile '$datafile' into table $tablename" );
	 print STDERR "imported $datafile\n";
       }

       $static_test_dbname = $dbname;
       $static_test_dbh = $db;

     }
   }
   return $static_test_dbh;
}


=head2 remove_test_db

  Args      : none
  Function  : drops test database which was created with test_db.
  Returntype: none
  Exceptions: none, will fail if no test database was created.
  Caller    : test scripts may remove the test database after usage,
              but should consider not, as upload takes some time.

=cut

sub remove_test_db {
  print STDERR "Dropping $static_test_dbname\n";
#  $static_test_dbh->do( "drop database $static_test_dbname" );
  $static_test_dbh->disconnect();
  $static_test_dbh = undef;
  $static_test_dbname = undef;
}


1;


__END__

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk
