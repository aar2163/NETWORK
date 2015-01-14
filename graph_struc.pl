#!/usr/bin/perl

use PDL;
use PDL::IO::Dumper;
#use PDL::Stats;
use lib '/home/andre/blib/lib/PDL/Lib/';
use lib '/home/andre/blib/arch/auto/PDL/Lib/Mylib/';
use Mylib;
use Time::HiRes qw(gettimeofday tv_interval);

use constant TRUE => 1;
use constant FALSE => 0;


$residue=-1;
@natom=();

$bEner=1;



if($#ARGV < 4)
{
 print "Usage: graph_struc.pl [i]input.pdb [i]adj.pld [i]enerd.pld [i]weight.pld [o]npath.pld\n";
 exit;
}

$adj = frestore($ARGV[1]);
$ener = frestore($ARGV[2]);
$ener = abs($ener);
#$ener2 = PDL::filter_ener($ener,$adj);
$ind1 = which($adj->flat == 0);
$ind2 = which($ener->flat > 1);
$ind3 = which($ener->flat != 0);
$intersect = intersect($ind1,$ind2);
$intersect2 = intersect($ind1,$ind3);

$nres = $adj->dice(0,X)->nelem;

#$d1 = $ener->flat->index($ind1);
$d1 = $ener->flat->index($intersect2);
$d2 = $ener->flat->index($intersect);

  ($mean,$prms,$median,$min,$max,$adev,$rms) = stats($d2);
#=begin
$d1 -= $mean;
$d1 /= (5*$prms);
$d1 += 1;
$d1 *= 0.5;

$d2 = $d1->where($d1 < 0.01);
$d2 .= 0;
$d2 = $d1->where($d1 > 0.99);
$d2 .= 0.99;
#=cut


#$ind2 = which($ener->flat < 40);
#$ind2 = which($ener2->flat == 0);
#$intersect = intersect($ind1,$ind2);

$weight = frestore($ARGV[3]);
#$com = frestore($ARGV[4]);

if($bEner)
{
 $weight .= 0.99;
 #$weight .= 0.5;
}

#$dice = $weight->flat->where(abs($weight->flat) < 1e-5);
#$dice .= 1e-5;

#$weight = -log(abs($weight));

#$weight = abs($weight);

#$max = $weight->max;

#$weight = $max-$weight+0.1;

#$dice = $weight->flat->index($intersect);
$d1 = $weight->flat->index($ind1);
$d2 = $ener->flat->index($ind1);

if($bEner)
{
 $d1 .= $d2;
}
else
{
 $d1 .= 0;
}


#$weight .= $adj;

 #$la = $weight->where($weight != 0);
 #print "$la\n";
 #exit;


#$diag = diagonal($adj,0,1);
#$diag .= 0;

#$adj /= $medres;


$dj = sumover($weight);
if(any $dj == 0)
{
 print "fatal error: $dj\n";
 exit;
}
#$diag = diagonal($dm,0,1);

#$diag .= $dj;

#$dinv = $dm->inv->sever;


$M = $weight / $dj->transpose;
$M = $M->transpose;


=begin
$H = zeroes($nres,$nres);
$H .= 1;
$h = diagonal($H,0,1);
$h .= 0;
=cut


$M = $M->sever;


=begin
 $time1 = [gettimeofday];
        $H2 = PDL::add($H,$M);
 $time2 = [gettimeofday];
  $elapsed = tv_interval($time1,$time2);
  print STDERR "$elapsed\n";


$Hr = sumover($H2);
#$Hr /= $Hr->nelem; ## NxN

#$Hr = zeroes($nres);
#PDL::cross_receive($H2,$Hr,253);

$Hr /= $Hr->maximum;
=cut

$npath = zeroes($nres);
#$dist = -log(abs($M)); #this is for conditional probability
$dist = $weight->copy;
$ind1 = which($weight->flat != 0);
$d1 = $dist->flat->index($ind1);
$d2 = $weight->flat->index($ind1);
$d1 .= -log(abs($d2));
#$dist = -log(abs($weight)); # this is for correlation coefficients

#$d1 = $weight->dice(292,X);
#print "$d1\n";
#exit;



=begin
$la1 = $weight->at(97,253+180);
$max = $weight->max;
$lala = whichND($weight == $max);
$dice = $weight->flat->where($weight->flat != 0);
$dice *= -1;
$dice += ($max+0.1);
#$dice = exp($dice);
$la2 = $weight->at(97,253+180);
print "$la1 $la2 $max $lala\n";
$dist = $weight;
exit;
=cut

      $npath = PDL::floyd($dist);
 $Hr = $npath/$npath->maximum;

fdump($npath,$ARGV[4]);

$dice = $weight->dice(97,X);
#$la2 = $ener->at(97,253+181);
$la1 = $npath->at(253+180);
$la2 = $npath->at(253+181);
$ind = which($dice != 0);
$la1 = $ener->dice(97,X);
$la2 = $la1->where($la1 > 0);
$ind = which($la1 > 0);
#print "$la2 $ind\n";
#exit;




$nfiles = $k;

$start{A} = 0;
$start{B} = $last{A};


$residue=-1;





    open (IN, "<$ARGV[0]");
    LINE: while (<IN>) {
       chomp;
       $line = $_;
       @data = split('',$_);
       $isatom = &join(0,3,@data);
       if($isatom eq  "ATOM")
       {
        $type = &join(13,15,@data);
        $type =~ s/^\s*(.*?)\s*$/$1/;
        $occup = &join(56,59,@data);
        if($type eq "N" or $type eq "C1")
        {
         $residue++;
        }
      
        $val = sclr($Hr->dice($residue));
        $val = sprintf("%3.2f",$val);

        @splited = split('\.', $occup);
 
         #### ************  ####
        $line =~ s/ $splited[0]\.$splited[1]/ $val/;
       }
         print "$line\n";
      }
    close IN;

sub join {
 my @array = @_;
 my $ii = shift(@array);
 my $jj = shift(@array);
 my $kk;
 my $string = $array[$ii];
 for ($kk=$ii+1;$kk<=$jj;$kk++)
 {
  $string = $string . "$array[$kk]";
 }
 return $string;
}

sub do_path {
  my $ii = $_[0];
  my $jj = $_[1];
  my $ii0 = $_[2];
  my $jj0 = $_[3];
  my $unw = $_[4];
  $nmat = (($unw == TRUE) ? $unmatrix : $nmatrix);
  $vmat = (($unw == TRUE) ? $uvmatrix : $vmatrix);
  my $kk = $nmat->at($ii,$jj);
   #$vm = $vmat->dice($ii);
   #$vm++;
  if($kk != -1)
  {
   @jjk = &get_index($jj);
   @iik = &get_index($ii);
   #printf "path $iik[0] $iik[1] $ii $jjk[0] $jjk[1] $jj kk: $kk ii0 $ii0 jj0 $jj0\n";
   do_path($kk,$jj,$ii0,$jj0,$unw);
   $vm = $vmat->dice($kk);
   #$vm = $vmat->dice($jj0);
   $vm++;
  }
  else
  {
   return;
  }
}


sub get_index {
 my $i = $_[0];
 $result[0] = $resindex[$i][0];
 $result[1] = $resindex[$i][1];
 return @result;
}


sub is_alpha
{
 my $ii = $_[0];
 my $ch = $_[1];
 my $jj,$jja,$jjb;

 my @listA = (110,119,153,163,87,94,58,72,32,40);
 my @listB = (14,26,184,199,67,76);
 @list = (($ch eq "A") ? @listA : @listB);
 for($jj=0;$jj<=($#list-1)/2;$jj++)
 {
  $jja = $list[2*$jj];
  $jjb = $list[(2*$jj)+1];
  if($ii >= $jja and $ii <= $jjb)
  {
   return 1;
  }
 }
 return 0;
}
sub is_beta
{
 my $ii = $_[0];
 my $ch = $_[1];
 my $jj,$jja,$jjb;

 my @listA = (167,172,198,201,5,12,46,51,76,80,137,141,126,134);
 my @listB = (2,6,45,48,79,83,172,175,166,169,135,141,122,126);
 @list = (($ch eq "A") ? @listA : @listB);
 for($jj=0;$jj<=($#list-1)/2;$jj++)
 {
  $jja = $list[2*$jj];
  $jjb = $list[(2*$jj)+1];
  if($ii >= $jja and $ii <= $jjb)
  {
   return 1;
  }
 }
 return 0;
}
sub is_pocket
{
 my $ii = $_[0];
 my $jj;

 my @pocket = (73,74,75,126,127,128,129,148,149,150);#,188,189,190,191,192,193,194,274,275,276,277,278,246,247,248,249);

 foreach $jj (@pocket)
 {
  if ($ii == $jj)
  {
   return 1;
  }
 }
 return 0;
}

sub which_group
{
 my $ii = $_[0];
 for($jj=0;$jj<$ngroup;$jj++)
 {
  if($ii >= $group[$jj][0] and $ii <= $group[$jj][1])
  {
   return $jj;
  }
 }
 return -1;
}
sub which_com
{
 my $num = $_[0];
 my $chain = $_[1];
 my $jj;
 for($jj=0;$jj<$ncom;$jj++)
 {
  foreach $elem (@{$community[$jj][$max_q]})
  {
   @iik = &get_index($elem);
   if($iik[0] == $num and $iik[1] eq $chain)
   {
    return $jj;
   }
  }
 }
 return -1;
}


sub split_blank {
#split entry
    my ($input) = shift;
    my @data;
    my (@array) = split (/ /, $input);
    my $k =0;
    foreach $var (@array) {
     if ($var =~ /\d/ or $var =~ /\w/ ) {
      $data[$k]=$var;
      $k++
     }
    }

    return @data;
}



