function basename(x) {
	gsub(".*/", "", x)
	return x
}

{
	split(basename($2), a, "_", seps)
	print a[1] "_" a[2] "\t" $1
}
