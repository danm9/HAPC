$(document).ready(function() {

// Loading samples 
// *******************************************************************************************************************

$('#load_sample').click(function() {
	hlas = 'B*51:01	2C\nB*51:01	3C\nB*35:01	7D\nB*58:01	5F';
	sequences = '001	A02:01:01G	A03:01:01G	B35:01:01G	B51:01:01G	C04:01:01G	C14:02:01G	ATGCGCTGCGGCTTCCGAGAC\n' +
				'002	A02:01:01G	A02:01:01G	B35:02:01G	B51:01:01G	C04:01:01G	C16:02:01G	ATGTGCTGCGGGATCCGAGAC';
	$('#hlas').text(hlas);
	$('#sequences').text(sequences);
	$('#hlas').val(hlas);
	$('#sequences').val(sequences);
});

// *******************************************************************************************************************

// Clearing
// **************************************************************

$('#clear').click(function() {
	$('#hlas').text('');
	$('#sequences').text('');
	$('#hlas').val('');
	$('#sequences').val('');
});

// **************************************************************

});
