cat << EOF  > "$archname"
#!/bin/sh
# This script was generated using Makeself $MS_VERSION


# MOD 0 START
# modifications for FESetup
# James Gebbie and Hannes H Loeffler (STFC Daresbury, UK) 2015
FES_install_top=\`pwd\`
FES_program='FESetup1.1'
FES_extract_dir='./'
# MOD 0 END

umask 077

CRCsum="$CRCsum"
MD5="$MD5sum"
TMPROOT=\${TMPDIR:=/tmp}
USER_PWD="\$PWD"; export USER_PWD

label="$LABEL"
script="$SCRIPT"
scriptargs="$SCRIPTARGS"
licensetxt="$LICENSE"
helpheader='$HELPHEADER'
targetdir="$archdirname"
filesizes="$filesizes"
keep="$KEEP"
quiet="n"

print_cmd_arg=""
if type printf > /dev/null; then
    print_cmd="printf"
elif test -x /usr/ucb/echo; then
    print_cmd="/usr/ucb/echo"
else
    print_cmd="echo"
fi

unset CDPATH

# MOD 1 START
# James Gebbie and Hannes H Loeffler
echo "Installing package \${label}..."
echo


FES_Control_C() {
    cd \$FES_install_top

    echo
    echo
    echo "Installation aborted by user. Cleaning up..."

    rm -rf \$targetdir

    exit 1
}

FES_error() {
    MS_Printf "\nERROR: \$1\n"
    MS_Printf "Aborting the installation of \$FES_program!\n"
    exit 1
}

FES_test_OS() {
    OS=\${1:-Linux}
    overwrite=\${2:-0}

    if type uname > /dev/null; then
        uname=uname
    else
	MS_Printf "WARNING: cannot determine OS type\n"
        uname=
    fi

    if [ x"\$overwrite" = x"1" -o -z "\$uname" ]; then
	MS_Printf "WARNING: system environent check skipped.\n"
	MS_Printf "Extracting all versions of \$FES_program.\n"
	FES_extract_dir='./'
	return
    fi

    MS_Printf "Determining system environment..."

    # What OS kernel are we running.
    kernel_name=\`uname -s\`

    if [ x"\$kernel_name" != x"\$OS" ]; then
	FES_error "unsupported operating system: \$kernel_name."
    fi

    # What Architecture are we running.	
    machine=\`uname -m\`

    case "\$machine" in
        i[3-6]86)
            FES_bits=32
            ;;
        x86_64|amd64)
            FES_bits=64
            ;;
        *)
            FES_error "unsupported architecture: \$machine."
    esac

    MS_Printf "    found \$kernel_name \$machine/\$FES_bits\n"

    # Extract the correct version.
    if [ "\$FES_bits" = "32" ]; then
	FES_extract_dir="./FESetup\$FES_bits"
    elif [ "\$FES_bits" = "64" ]; then
	FES_extract_dir="./FESetup\$FES_bits"
    fi
}

FES_test_python() {
    FES_python=\${1:-python}
    python_required=\${2:-2.7}
    overwrite=\${3:-0}

    if [ x"\$overwrite" = x"1" ]; then
	MS_Printf "WARNING: Python version check skipped.\n"
	MS_Printf "Extracting all versions of \$FES_program.\n"
	return
    fi

    MS_Printf "Checking for Python version \$python_required..."

    python_version=\`\$FES_python --version 2>&1\`

    if [ \$? -ne 0 ]; then
        FES_error "Python interpreter \$FES_python cannot be run."
    fi

    python_version=\`echo \$python_version | cut -c8-10\`

    if [ "\$python_version" != "\$python_required" ]; then
        FES_error "Incompatible version \$python_version found."
    else 
        MS_Printf "   found \$python_version\n"
    fi
}

FES_make_script() {
    scr=FESetup
    tmpl=FESetup.tmpl
    curr_dir=\`pwd -P\`

    if ! cd FESetup\$1/bin > /dev/null; then
	FES_error "directory FESetup\$1/bin does not exist."
    fi

    if [ ! -f \$tmpl ]; then
        FES_error "template file \$tmpl not found."
    fi

    sed -e "s,%TOPDIR%,\$FES_TOP/FESetup\$1," \\
	-e "s,%PYINP%,\$FES_python," \\
	-e "s,%AMBER%,\$FES_AMBER," \\
	-e "s,%GROMACS%,\$FES_GROMACS," \\
	-e "s,%NAMD%,\$FES_NAMD," \\
	-e "s,%DLPOLY%,\$FES_DLPOLY," \$tmpl > \$scr

    chmod +x \$scr
    rm -f \$tmpl

    rm -f dGprep.py
    ln -sf ../lib/python2.7/site-packages/FESetup/ui/dGprep.py

    cd \$curr_dir
}

FES_check_path () {
    path=\$1
    exes=\$2

    if [ ! -d "\$path" -o -z "\$path" ]; then
        FES_check_path_ok=0
        return
    fi

    found=0
    total=0

    for exe in \$exes; do
        exes2=\$path/\$exe

        for exe2 in \$exes2; do
            total=`expr \$total + 1`

            if [ -x "\$exe2" ]; then
                found=`expr \$found + 1`
            fi 
        done
    done

    if [ \$found -ne \$total ]; then
        FES_check_path_ok=0
    else
        FES_check_path_ok=1
    fi
}

FES_post_install() {
    FES_TOP=\`pwd -P\`

    if [ -z "\$FES_bits" ]; then
        bits="32 or 64"
    else
        bits=\$FES_bits
    fi 

    cat <<_FES_POST

=========================  post installation setup =========================

The installer will now ask for the names of several directories where the
MD software packages AMBER, GROMACS, NAMD and DL_POLY are installed.
All options either have defaults or are optional.  You can change this
any time later in \$targetdir/FESetup\$bits/bin/FESetup.

Confirm defaults by pressing the Enter key.

This installer comes with the AMBER tools required for \$FES_program but you
can choose to point AMBERHOME to your own installation.  If you choose your
own please be aware that you will need to modify dat/antechamber/PARMCHK.DAT
as contained in this installer package to avoid segmentation faults with
parmchk 2.

_FES_POST

    amber=n
    amber_exe_list="sander antechamber teLeap sqm"

    if [ -n "\$AMBERHOME" ]; then
	MS_Printf "\\\$AMBERHOME is set to \$AMBERHOME, do you wish to use this? [y/N] "

	read amber
    fi

    case "\$amber" in
	y*|Y*)
            FES_AMBER=\$AMBERHOME

	    FES_check_path \$FES_AMBER/bin "\$amber_exe_list"

	    if [ \$FES_check_path_ok -lt 1 ]; then
		echo 'WARNING: AMBER not properly installed, using default path!'
		FES_AMBER='\$FES_HOME/amber14'
	    fi
	    ;;
	*)
            MS_Printf 'AMBER directory? [\$FES_HOME/amber14] '

	    read amber2

	    if [ -z "\$amber2" ]; then
		FES_AMBER='\$FES_HOME/amber14'
	    else
		FES_AMBER=\$amber2

		FES_check_path \$FES_AMBER/bin "\$amber_exe_list"

		if [ \$FES_check_path_ok -lt 1 ]; then
		    echo 'WARNING: AMBER not properly installed, '\
                         'using default path!'
		    FES_AMBER='\$FES_HOME/amber14'
		fi
	    fi
    esac


    gromacs=n

    if [ -n "\$GMXBIN" ]; then
	gmxbin=\`dirname \$GMXBIN\`
	
	MS_Printf "GROMACS appears to be installed in \$gmxbin, do you wish to use this? [Y/n] "

	read gromacs
    fi

    case "\$gromacs" in
	n*|N*)
            MS_Printf 'GROMACS directory? [] '

	    read gromacs2

	    if [ -n "\$gromacs2" ]; then
		FES_GROMACS=\$gromacs2
	    fi
	    ;;
	*)
            FES_GROMACS=\$gmxbin
    esac

    if [ -n "\$FES_GROMACS" ]; then
	FES_check_path \$FES_GROMACS/bin "mdrun* grompp*"
	
	if [ \$FES_check_path_ok -lt 1 ]; then
	    echo 'WARNING: GROMACS not properly installed, using empty path!'
	    FES_GROMACS=
	fi
    fi


    MS_Printf 'NAMD directory (path contains namd2)? [] '

    read namd

    if [ -n "\$namd" ]; then
	FES_NAMD=\$namd
	
	FES_check_path \$FES_NAMD "namd2"

	if [ \$FES_check_path_ok -lt 1 ]; then
	    echo 'WARNING: NAMD not properly installed, using empty path!'
	    FES_NAMD=
	fi
    fi


    MS_Printf 'DL_POLY directory (path contains execute/DLPOLY.Z)? [] '

    read dlpoly

    if [ -n "\$dlpoly" ]; then
	FES_DLPOLY=\$dlpoly

	FES_check_path \$FES_DLPOLY/execute "DLPOLY.Z"

	if [ \$FES_check_path_ok -lt 1 ]; then
	    echo 'WARNING: DL_POLY not properly installed, using empty path!'
	    FES_DLPOLY=
	fi
    fi


    if [ -z "\$FES_bits" ]; then
	FES_make_script 32
	FES_make_script 64
    else
	FES_make_script \$FES_bits
    fi

    if [ "\$FES_bits" -eq 64 -o -z "\$FES_bits" ]; then
	ucs=\`\$FES_python -c 'import sys; print sys.maxunicode'\`

	if [ "\$ucs" -gt 65535 ]; then
	    nbytes=4
	else
	    nbytes=2
	fi

	curr_dir=\`pwd -P\`

	if ! cd FESetup64/lib > /dev/null; then
	    FES_error "directory FESetup64/lib does not exist."
	fi

	if [ ! -f libboost_python.so.1.57.0.UCS\$nbytes ]; then
	    FES_error "libboost_python.so.1.57.0.UCS\$nbytes does not exist."
	fi

	rm -f libboost_python.so.1.57.0
	ln -sf libboost_python.so.1.57.0.UCS\$nbytes libboost_python.so.1.57.0

	if ! cd python2.7/site-packages > /dev/null; then
	    FES_error "directory python2.7/site-packages does not exist."
	fi

	if [ ! -d Sire.UCS\$nbytes ]; then
	    FES_error "direcotry Sire.UCS\$nbytes does not exist."
	fi

	if [ ! -d numpy.UCS\$nbytes ]; then
	    FES_error "direcotry numpy.UCS\$nbytes does not exist."
	fi

	rm -f Sire numpy
	ln -sf Sire.UCS\$nbytes Sire
	ln -sf numpy.UCS\$nbytes numpy

	cd \$curr_dir
    fi
}
# MOD 1 END

MS_Printf()
{
    \$print_cmd \$print_cmd_arg "\$1"
}

MS_PrintLicense()
{
  if test x"\$licensetxt" != x; then
    echo "\$licensetxt"
    while true
    do
      MS_Printf "Please type y to accept, n otherwise: "
      read yn
      if test x"\$yn" = xn; then
        keep=n
	eval \$finish; exit 1
        break;
      elif test x"\$yn" = xy; then
        break;
      fi
    done
  fi
}

MS_diskspace()
{
	(
	if test -d /usr/xpg4/bin; then
		PATH=/usr/xpg4/bin:\$PATH
	fi
	df -kP "\$1" | tail -1 | awk '{ if (\$4 ~ /%/) {print \$3} else {print \$4} }'
	)
}

MS_dd()
{
    blocks=\`expr \$3 / 1024\`
    bytes=\`expr \$3 % 1024\`
    dd if="\$1" ibs=\$2 skip=1 obs=1024 conv=sync 2> /dev/null | \\
    { test \$blocks -gt 0 && dd ibs=1024 obs=1024 count=\$blocks ; \\
      test \$bytes  -gt 0 && dd ibs=1 obs=1024 count=\$bytes ; } 2> /dev/null
}

MS_dd_Progress()
{
    if test x"\$noprogress" = xy; then
        MS_dd \$@
        return \$?
    fi
    file="\$1"
    offset=\$2
    length=\$3
    pos=0
    bsize=4194304
    while test \$bsize -gt \$length; do
        bsize=\`expr \$bsize / 4\`
    done
    blocks=\`expr \$length / \$bsize\`
    bytes=\`expr \$length % \$bsize\`
    (
        dd ibs=\$offset skip=1 2>/dev/null
        pos=\`expr \$pos \+ \$bsize\`
        MS_Printf "     0%% " 1>&2
        if test \$blocks -gt 0; then
            while test \$pos -le \$length; do
                dd bs=\$bsize count=1 2>/dev/null
                pcent=\`expr \$length / 100\`
                pcent=\`expr \$pos / \$pcent\`
                if test \$pcent -lt 100; then
                    MS_Printf "\b\b\b\b\b\b\b" 1>&2
                    if test \$pcent -lt 10; then
                        MS_Printf "    \$pcent%% " 1>&2
                    else
                        MS_Printf "   \$pcent%% " 1>&2
                    fi
                fi
                pos=\`expr \$pos \+ \$bsize\`
            done
        fi
        if test \$bytes -gt 0; then
            dd bs=\$bytes count=1 2>/dev/null
        fi
        MS_Printf "\b\b\b\b\b\b\b" 1>&2
        MS_Printf " 100%%  " 1>&2
    ) < "\$file"
}

MS_Help()
{
    cat << EOH >&2
\${helpheader}Makeself version $MS_VERSION
 1) Getting help or info about \$0 :
  \$0 --help   Print this message
  \$0 --info   Print embedded info : title, default target directory, embedded script ...
  \$0 --lsm    Print embedded lsm entry (or no LSM)
  \$0 --list   Print the list of files in the archive
  \$0 --check  Checks integrity of the archive

 2) Running \$0 :
  \$0 [options] [--] [additional arguments to embedded script]
  with following options (in that order)
  --confirm             Ask before running embedded script
  --quiet		Do not print anything except error messages
  --noexec              Do not run embedded script
  --keep                Do not erase target directory after running
			the embedded script
  --noprogress          Do not show the progress during the decompression
  --nox11               Do not spawn an xterm
  --nochown             Do not give the extracted files to the current user
  --target dir          Extract directly to a target directory
                        directory path can be either absolute or relative
  --tar arg1 [arg2 ...] Access the contents of the archive through the tar command
  --extract-all         FESetup: extract complete content, skip all checks
  --python interpreter  FESetup: path to alternative Python interpreter
  --                    Following arguments will be passed to the embedded script
EOH
}

MS_Check()
{
    OLD_PATH="\$PATH"
    PATH=\${GUESS_MD5_PATH:-"\$OLD_PATH:/bin:/usr/bin:/sbin:/usr/local/ssl/bin:/usr/local/bin:/opt/openssl/bin"}
	MD5_ARG=""
    MD5_PATH=\`exec <&- 2>&-; which md5sum || type md5sum\`
    test -x "\$MD5_PATH" || MD5_PATH=\`exec <&- 2>&-; which md5 || type md5\`
	test -x "\$MD5_PATH" || MD5_PATH=\`exec <&- 2>&-; which digest || type digest\`
    PATH="\$OLD_PATH"

    if test x"\$quiet" = xn; then
		MS_Printf "Verifying archive integrity..."
    fi
    offset=\`head -n $SKIP "\$1" | wc -c | tr -d " "\`
    verb=\$2
    i=1
    for s in \$filesizes
    do
		crc=\`echo \$CRCsum | cut -d" " -f\$i\`
		if test -x "\$MD5_PATH"; then
			if test x"\`basename \$MD5_PATH\`" = xdigest; then
				MD5_ARG="-a md5"
			fi
			md5=\`echo \$MD5 | cut -d" " -f\$i\`
			if test x"\$md5" = x00000000000000000000000000000000; then
				test x"\$verb" = xy && echo " \$1 does not contain an embedded MD5 checksum." >&2
			else
				md5sum=\`MS_dd "\$1" \$offset \$s | eval "\$MD5_PATH \$MD5_ARG" | cut -b-32\`;
				if test x"\$md5sum" != x"\$md5"; then
					echo "Error in MD5 checksums: \$md5sum is different from \$md5" >&2
					exit 2
				else
					test x"\$verb" = xy && MS_Printf " MD5 checksums are OK." >&2
				fi
				crc="0000000000"; verb=n
			fi
		fi
		if test x"\$crc" = x0000000000; then
			test x"\$verb" = xy && echo " \$1 does not contain a CRC checksum." >&2
		else
			sum1=\`MS_dd "\$1" \$offset \$s | CMD_ENV=xpg4 cksum | awk '{print \$1}'\`
			if test x"\$sum1" = x"\$crc"; then
				test x"\$verb" = xy && MS_Printf " CRC checksums are OK." >&2
			else
				echo "Error in checksums: \$sum1 is different from \$crc" >&2
				exit 2;
			fi
		fi
		i=\`expr \$i + 1\`
		offset=\`expr \$offset + \$s\`
    done
    if test x"\$quiet" = xn; then
		echo " All good."
    fi
}

UnTAR()
{
    if test x"\$quiet" = xn; then
# MOD 2 START
# James Gebbie and Hannes H Loeffler
	tar \$1vf - \$2 2>&1 || { echo Extraction failed. > /dev/tty; kill -15 \$$; }
    else
	tar \$1f - \$2 2>&1 || { echo Extraction failed. > /dev/tty; kill -15 \$$;
}
# MOD 2 END
    fi
}

finish=true
xterm_loop=
noprogress=$NOPROGRESS
nox11=$NOX11
copy=$COPY
ownership=y
verbose=n

initargs="\$@"

while true
do
    case "\$1" in
    -h | --help)
	MS_Help
	exit 0
	;;
    -q | --quiet)
	quiet=y
	noprogress=y
	shift
	;;
    --info)
	echo Identification: "\$label"
	echo Target directory: "\$targetdir"
	echo Uncompressed size: $USIZE KB
	echo Compression: $COMPRESS
	echo Date of packaging: $DATE
	echo Built with Makeself version $MS_VERSION on $OSTYPE
	echo Build command was: "$MS_COMMAND"
	if test x"\$script" != x; then
	    echo Script run after extraction:
	    echo "    " \$script \$scriptargs
	fi
	if test x"$copy" = xcopy; then
		echo "Archive will copy itself to a temporary location"
	fi
	if test x"$KEEP" = xy; then
	    echo "directory \$targetdir is permanent"
	else
	    echo "\$targetdir will be removed after extraction"
	fi
	exit 0
	;;
    --dumpconf)
	echo LABEL=\"\$label\"
	echo SCRIPT=\"\$script\"
	echo SCRIPTARGS=\"\$scriptargs\"
	echo archdirname=\"$archdirname\"
	echo KEEP=$KEEP
	echo COMPRESS=$COMPRESS
	echo filesizes=\"\$filesizes\"
	echo CRCsum=\"\$CRCsum\"
	echo MD5sum=\"\$MD5\"
	echo OLDUSIZE=$USIZE
	echo OLDSKIP=`expr $SKIP + 1`
	exit 0
	;;
    --lsm)
cat << EOLSM
EOF
eval "$LSM_CMD"
cat << EOF  >> "$archname"
EOLSM
	exit 0
	;;
    --list)
	echo Target directory: \$targetdir
	offset=\`head -n $SKIP "\$0" | wc -c | tr -d " "\`
	for s in \$filesizes
	do
	    MS_dd "\$0" \$offset \$s | eval "$GUNZIP_CMD" | UnTAR t
	    offset=\`expr \$offset + \$s\`
	done
	exit 0
	;;
	--tar)
	offset=\`head -n $SKIP "\$0" | wc -c | tr -d " "\`
	arg1="\$2"
    if ! shift 2; then MS_Help; exit 1; fi
	for s in \$filesizes
	do
	    MS_dd "\$0" \$offset \$s | eval "$GUNZIP_CMD" | tar "\$arg1" - "\$@"
	    offset=\`expr \$offset + \$s\`
	done
	exit 0
	;;
    --check)
	MS_Check "\$0" y
	exit 0
	;;
    --confirm)
	verbose=y
	shift
	;;
	--noexec)
	script=""
	shift
	;;
    --keep)
	keep=y
	shift
	;;
    --target)
	keep=y
	targetdir=\${2:-.}
    if ! shift 2; then MS_Help; exit 1; fi
	;;
    --noprogress)
	noprogress=y
	shift
	;;
    --nox11)
	nox11=y
	shift
	;;
    --nochown)
	ownership=n
	shift
	;;
    --xwin)
	if test "$NOWAIT" = n; then
		finish="echo Press Return to close this window...; read junk"
	fi
	xterm_loop=1
	shift
	;;
    --phase2)
	copy=phase2
	shift
	;;
# MOD 3 START
# James Gebbie and Hannes H Loeffler
    --extract-all)
        extract_all=1
	shift
	;;
    --python)
        python_interpreter=\${2:-.}
	if ! shift 2; then MS_Help; exit 1; fi
        ;;
# MOD 3 END
    --)
	shift
	break ;;
    -*)
	echo Unrecognized flag : "\$1" >&2
	MS_Help
	exit 1
	;;
    *)
	break ;;
    esac
done

# MOD 4 START
# James Gebbie and Hannes H Loeffler
trap FES_Control_C 2
FES_test_OS 'Linux' \$extract_all
FES_test_python "\$python_interpreter" '2.7' \$extract_all
# MOD 4 END

if test x"\$quiet" = xy -a x"\$verbose" = xy; then
	echo Cannot be verbose and quiet at the same time. >&2
	exit 1
fi

if test x"\$copy" \!= xphase2; then
    MS_PrintLicense
fi

case "\$copy" in
copy)
    tmpdir=\$TMPROOT/makeself.\$RANDOM.\`date +"%y%m%d%H%M%S"\`.\$\$
    mkdir "\$tmpdir" || {
	echo "Could not create temporary directory \$tmpdir" >&2
	exit 1
    }
    SCRIPT_COPY="\$tmpdir/makeself"
    echo "Copying to a temporary location..." >&2
    cp "\$0" "\$SCRIPT_COPY"
    chmod +x "\$SCRIPT_COPY"
    cd "\$TMPROOT"
    exec "\$SCRIPT_COPY" --phase2 -- \$initargs
    ;;
phase2)
    finish="\$finish ; rm -rf \`dirname \$0\`"
    ;;
esac

if test x"\$nox11" = xn; then
    if tty -s; then                 # Do we have a terminal?
	:
    else
        if test x"\$DISPLAY" != x -a x"\$xterm_loop" = x; then  # No, but do we have X?
            if xset q > /dev/null 2>&1; then # Check for valid DISPLAY variable
                GUESS_XTERMS="xterm gnome-terminal rxvt dtterm eterm Eterm xfce4-terminal lxterminal kvt konsole aterm terminology"
                for a in \$GUESS_XTERMS; do
                    if type \$a >/dev/null 2>&1; then
                        XTERM=\$a
                        break
                    fi
                done
                chmod a+x \$0 || echo Please add execution rights on \$0
                if test \`echo "\$0" | cut -c1\` = "/"; then # Spawn a terminal!
                    exec \$XTERM -title "\$label" -e "\$0" --xwin "\$initargs"
                else
                    exec \$XTERM -title "\$label" -e "./\$0" --xwin "\$initargs"
                fi
            fi
        fi
    fi
fi

if test x"\$targetdir" = x.; then
    tmpdir="."
else
    if test x"\$keep" = xy; then
	if test x"\$quiet" = xn; then
	    echo "Creating directory \$targetdir" >&2
	fi
	tmpdir="\$targetdir"
	dashp="-p"
    else
	tmpdir="\$TMPROOT/selfgz\$\$\$RANDOM"
	dashp=""
    fi

    mkdir \$dashp \$tmpdir || {
	echo 'Cannot create target directory' \$tmpdir >&2
	echo 'You should try option --target dir' >&2
	eval \$finish
	exit 1
    }
fi

location="\`pwd\`"
if test x"\$SETUP_NOCHECK" != x1; then
    MS_Check "\$0"
fi
offset=\`head -n $SKIP "\$0" | wc -c | tr -d " "\`

if test x"\$verbose" = xy; then
	MS_Printf "About to extract $USIZE KB in \$tmpdir ... Proceed ? [Y/n] "
	read yn
	if test x"\$yn" = xn; then
		eval \$finish; exit 1
	fi
fi

if test x"\$quiet" = xn; then
	MS_Printf "Uncompressing \$label"
fi
res=3
if test x"\$keep" = xn; then
    trap 'echo Signal caught, cleaning up >&2; cd \$TMPROOT; /bin/rm -rf \$tmpdir; eval \$finish; exit 15' 1 2 3 15
fi

leftspace=\`MS_diskspace \$tmpdir\`
if test -n "\$leftspace"; then
    if test "\$leftspace" -lt $USIZE; then
        echo
        echo "Not enough space left in "\`dirname \$tmpdir\`" (\$leftspace KB) to decompress \$0 ($USIZE KB)" >&2
        if test x"\$keep" = xn; then
            echo "Consider setting TMPDIR to a directory with more free space."
        fi
        eval \$finish; exit 1
    fi
fi

for s in \$filesizes
do
# MOD 5 START
# James Gebbie and Hannes H Loeffler
    if MS_dd_Progress "\$0" \$offset \$s | eval "$GUNZIP_CMD" | ( cd "\$tmpdir"; UnTAR xp \$FES_extract_dir) 1>/dev/null; then
# MOD 5 END
		if test x"\$ownership" = xy; then
			(PATH=/usr/xpg4/bin:\$PATH; cd "\$tmpdir"; chown -R \`id -u\` .;  chgrp -R \`id -g\` .)
		fi
    else
		echo >&2
		echo "Unable to decompress \$0" >&2
		eval \$finish; exit 1
    fi
    offset=\`expr \$offset + \$s\`
done
if test x"\$quiet" = xn; then
	echo
fi

cd "\$tmpdir"
res=0
if test x"\$script" != x; then
    if test x"\$verbose" = x"y"; then
		MS_Printf "OK to execute: \$script \$scriptargs \$* ? [Y/n] "
		read yn
		if test x"\$yn" = x -o x"\$yn" = xy -o x"\$yn" = xY; then
			eval "\"\$script\" \$scriptargs \"\\\$@\""; res=\$?;
		fi
    else
		eval "\"\$script\" \$scriptargs \"\\\$@\""; res=\$?
    fi
    if test "\$res" -ne 0; then
		test x"\$verbose" = xy && echo "The program '\$script' returned an error code (\$res)" >&2
    fi
fi
if test x"\$keep" = xn; then
    cd \$TMPROOT
    /bin/rm -rf \$tmpdir
fi

# MOD 6 START
# James Gebbie and Hannes H Loeffler
FES_post_install

cat <<_FES_POST

========== Installation of \$label is now complete. ==========

_FES_POST
# MOD 6 END

eval \$finish; exit \$res
EOF
